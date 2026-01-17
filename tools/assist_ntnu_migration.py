import os
import re
import csv
import json
import numpy as np

# Paths
SOURCE_BASE = "/home/shug8104/alm/Actuator-Line-Code/tutorials/NTNU_ALM/constant"
DEST_BASE = "/home/shug8104/FVW-CPP/data/NTNU"
TURBINE_PROPS = os.path.join(SOURCE_BASE, "turbineProperties", "T2")
AIRFOIL_DIR = os.path.join(SOURCE_BASE, "airfoilProperties")

def ensure_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def parse_openfoam_list(content, list_name):
    # Find the list by name, handling newlines and whitespace
    pattern = list_name + r"\s*\n?\s*\(\s*([\s\S]*?)\s*\)\s*;"
    match = re.search(pattern, content)
    if not match:
        print(f"Warning: Could not find list {list_name}")
        return []
    
    inner_content = match.group(1)
    # Remove C++ style comments
    inner_content = re.sub(r'//.*', '', inner_content)
    
    # Check if it's a list of strings (like Airfoils) or list of tuples
    if '"' in inner_content:
        # List of strings
        items = re.findall(r'"([^"]+)"', inner_content)
        return items
    else:
        # List of tuples/numbers
        # Remove parens and extra spaces
        items = []
        # Split by matching sets of parentheses for rows
        rows = re.findall(r'\(([^)]+)\)', inner_content)
        for row in rows:
            vals = [float(x) for x in row.split()]
            items.append(vals)
        return items

def read_turbine_props():
    with open(TURBINE_PROPS, 'r') as f:
        content = f.read()
    
    props = {}
    
    # Simple scalar values
    scalar_map = {
        'NumBl': ('nBlades', int),
        'TipRad': ('rTip', float),
        'HubRad': ('rHub', float),
        'TowerHt': ('hubHeight', float)
    }
    
    for foam_key, (json_key, dtype) in scalar_map.items():
        match = re.search(rf"{foam_key}\s+([0-9.]+);", content)
        if match:
            props[json_key] = dtype(match.group(1))
    
    # Lists
    props['Airfoils'] = parse_openfoam_list(content, 'Airfoils')
    props['BladeData'] = parse_openfoam_list(content, 'BladeData')
    props['rotorDesign'] = parse_openfoam_list(content, 'rotorDesign')
    
    return props

def write_turbine_params(props):
    data = {
        "rTip": props.get('rTip', 0.447),
        "rHub": props.get('rHub', 0.045),
        "nBlades": props.get('nBlades', 3),
        "hubHeight": props.get('hubHeight', 0.8)
    }
    
    with open(os.path.join(DEST_BASE, "turbine_params.json"), 'w') as f:
        json.dump(data, f, indent=4)
    print(f"Created turbine_params.json")

def write_blade_geometry(props):
    # rotorDesign: radius(m), c(m), twist(deg), thickness-ratio
    # BladeData: radius(m), Re, Tu, airfoil_index
    
    # We need to merge these. 
    # The rotorDesign usually has the geometry.
    # The BladeData maps radius to airfoil index.
    # We will interpolate airfoil index into rotorDesign radial stations if needed, 
    # or just use rotorDesign stations and find nearest airfoil index.
    
    rotor_design = props['rotorDesign'] # [r, c, twist, t/c]
    blade_data = props['BladeData']     # [r, Re, Tu, airfoil_idx]
    
    # Create interpolator for airfoil index
    bd_r = [row[0] for row in blade_data]
    bd_idx = [row[3] for row in blade_data]
    
    with open(os.path.join(DEST_BASE, "blade_geometry.csv"), 'w') as f:
        writer = csv.writer(f)
        writer.writerow(["r", "chord", "twist", "airfoil_index"])
        
        for row in rotor_design:
            r, c, twist, _ = row
            
            # Find nearest airfoil index (step function usually)
            # Simple nearest neighbor for now or interp? 
            # Airfoil indices are integers, so usually it's defined by zones.
            # Let's find the closest r in blade_data and use that index.
            
            # Find index in bd_r closest to r
            idx = int(bd_idx[np.abs(np.array(bd_r) - r).argmin()])
            
            writer.writerow([r, c, twist, idx])
            
    print(f"Created blade_geometry.csv")

def write_airfoil_list(props):
    airfoils = props['Airfoils']
    with open(os.path.join(DEST_BASE, "airfoil_list.txt"), 'w') as f:
        for af in airfoils:
            f.write(f"{af}\n")
    print(f"Created airfoil_list.txt")
    return airfoils

def convert_airfoils(airfoil_names):
    dest_dir = os.path.join(DEST_BASE, "airfoils")
    ensure_dir(dest_dir)
    
    for name in airfoil_names:
        src_path = os.path.join(AIRFOIL_DIR, name)
        dest_path = os.path.join(dest_dir, name + ".dat")
        
        if not os.path.exists(src_path):
            print(f"Warning: Source airfoil {src_path} not found")
            continue
            
        with open(src_path, 'r') as f:
            content = f.read()
            
        # Parse OpenFOAM table
        # airfoilData ( ... )
        data = parse_openfoam_list(content, 'airfoilData')
        
        with open(dest_path, 'w') as f:
            f.write(f"Airfoil {name} converted from OpenFOAM\n")
            f.write("Generated by migration script\n")
            f.write("   1        Number of airfoil tables\n")
            f.write("   0.0      Table ID\n")
            f.write("   0.0      Stall angle\n")
            f.write("   0.0      Unused\n")
            f.write("   0.0      Unused\n")
            f.write("   0.0      Unused\n")
            f.write("   0.0      Zero Cn angle\n")
            f.write("   0.0      Cn slope\n")
            f.write("   0.0      Cn pos stall\n")
            f.write("   0.0      Cn neg stall\n")
            f.write("   0.0      Min Cd ang\n")
            f.write("   0.0      Min Cd val\n")
            
            for row in data:
                alpha = row[0]
                cl = row[1]
                cd = row[2]
                cm = 0.0 # Missing in source
                f.write(f"{alpha:8.3f} {cl:8.4f} {cd:8.4f} {cm:8.4f}\n")
                
    print(f"Converted {len(airfoil_names)} airfoils")

def main():
    print("Starting NTNU migration...")
    ensure_dir(DEST_BASE)
    
    props = read_turbine_props()
    # print("Props:", props.keys())
    
    write_turbine_params(props)
    write_blade_geometry(props)
    airfoils = write_airfoil_list(props)
    convert_airfoils(airfoils)
    
    print("Migration complete!")

if __name__ == "__main__":
    main()
