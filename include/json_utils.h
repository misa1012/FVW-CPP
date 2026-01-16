#ifndef JSON_UTILS_H
#define JSON_UTILS_H

#include <string>
#include <vector>
#include <map>
#include <variant>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <stdexcept>

// A minimal JSON parser for FVW-CPP configuration
// Supports: Object, Array, String, Number, Boolean, Null
// Limitation: No Unicode escape support, basic error handling

namespace fvw {
namespace json {

enum class GenericType {
    Null,
    Boolean,
    Number,
    String,
    Array,
    Object
};

struct Value;

using NullType = std::monostate;
using BooleanType = bool;
using NumberType = double;
using StringType = std::string;
using ArrayType = std::vector<Value>;
using ObjectType = std::map<std::string, Value>;

struct Value {
    GenericType type;
    std::variant<NullType, BooleanType, NumberType, StringType, ArrayType, ObjectType> data;

    Value() : type(GenericType::Null), data(NullType{}) {}
    Value(bool b) : type(GenericType::Boolean), data(b) {}
    Value(double d) : type(GenericType::Number), data(d) {}
    Value(int i) : type(GenericType::Number), data(static_cast<double>(i)) {}
    Value(const std::string& s) : type(GenericType::String), data(s) {}
    Value(const char* s) : type(GenericType::String), data(std::string(s)) {}
    Value(const ArrayType& a) : type(GenericType::Array), data(a) {}
    Value(const ObjectType& o) : type(GenericType::Object), data(o) {}

    // Helpers
    bool as_bool() const {
        if (type != GenericType::Boolean) throw std::runtime_error("Type mismatch: expected boolean");
        return std::get<BooleanType>(data);
    }
    double as_double() const {
        if (type != GenericType::Number) throw std::runtime_error("Type mismatch: expected number");
        return std::get<NumberType>(data);
    }
    int as_int() const {
        return static_cast<int>(as_double());
    }
    std::string as_string() const {
        if (type != GenericType::String) throw std::runtime_error("Type mismatch: expected string");
        return std::get<StringType>(data);
    }
    const ArrayType& as_array() const {
        if (type != GenericType::Array) throw std::runtime_error("Type mismatch: expected array");
        return std::get<ArrayType>(data);
    }
    const ObjectType& as_object() const {
        if (type != GenericType::Object) throw std::runtime_error("Type mismatch: expected object");
        return std::get<ObjectType>(data);
    }
    bool contains(const std::string& key) const {
        if(type != GenericType::Object) return false;
        const auto& obj = std::get<ObjectType>(data);
        return obj.find(key) != obj.end();
    }
    const Value& operator[](const std::string& key) const {
        const auto& obj = as_object();
        auto it = obj.find(key);
        if (it == obj.end()) throw std::runtime_error("Key not found: " + key);
        return it->second;
    }
};

class Parser {
    std::string_view src;
    size_t pos = 0;

    char peek() {
        skip_whitespace();
        if (pos >= src.size()) return '\0';
        return src[pos];
    }

    char advance() {
        char c = peek();
        if (pos < src.size()) pos++;
        return c;
    }

    void skip_whitespace() {
        while (pos < src.size() && std::isspace(src[pos])) pos++;
    }

    void expect(char c) {
        if (peek() != c) throw std::runtime_error(std::string("Expected '") + c + "'");
        pos++;
    }

    std::string parse_string_raw() {
        expect('"');
        std::string res;
        while (pos < src.size()) {
            char c = src[pos++];
            if (c == '"') return res;
            if (c == '\\') {
                if (pos >= src.size()) throw std::runtime_error("Unexpected end in string escape");
                char esc = src[pos++];
                if (esc == '"') res += '"';
                else if (esc == '\\') res += '\\';
                else if (esc == '/') res += '/';
                else if (esc == 'b') res += '\b';
                else if (esc == 'f') res += '\f';
                else if (esc == 'n') res += '\n';
                else if (esc == 'r') res += '\r';
                else if (esc == 't') res += '\t';
                // Unicode \uXXXX not implemented for simplicity
                else res += esc;
            } else {
                res += c;
            }
        }
        throw std::runtime_error("Unterminated string");
    }

    Value parse_value() {
        char c = peek();
        if (c == 'n') {
            if (src.substr(pos, 4) == "null") { pos += 4; return Value(); }
        }
        if (c == 't') {
            if (src.substr(pos, 4) == "true") { pos += 4; return Value(true); }
        }
        if (c == 'f') {
            if (src.substr(pos, 5) == "false") { pos += 5; return Value(false); }
        }
        if (c == '"') {
            return Value(parse_string_raw());
        }
        if (c == '[') {
            expect('[');
            ArrayType arr;
            if (peek() != ']') {
                while (true) {
                    arr.push_back(parse_value());
                    if (peek() == ',') advance();
                    else break;
                }
            }
            expect(']');
            return Value(arr);
        }
        if (c == '{') {
            expect('{');
            ObjectType obj;
            if (peek() != '}') {
                while (true) {
                    std::string key = parse_string_raw();
                    expect(':');
                    obj[key] = parse_value();
                    if (peek() == ',') advance();
                    else break;
                }
            }
            expect('}');
            return Value(obj);
        }
        if (c == '-' || std::isdigit(c)) {
            size_t start = pos;
            if (c == '-') pos++;
            while (pos < src.size() && std::isdigit(src[pos])) pos++;
            if (pos < src.size() && src[pos] == '.') {
                pos++;
                while (pos < src.size() && std::isdigit(src[pos])) pos++;
            }
            if (pos < src.size() && (src[pos] == 'e' || src[pos] == 'E')) {
                pos++;
                if (pos < src.size() && (src[pos] == '+' || src[pos] == '-')) pos++;
                while (pos < src.size() && std::isdigit(src[pos])) pos++;
            }
            std::string num_str(src.substr(start, pos - start));
            return Value(std::stod(num_str));
        }
        throw std::runtime_error(std::string("Unexpected character: ") + c);
    }

public:
    Value parse(std::string_view source) {
        src = source;
        pos = 0;
        return parse_value();
    }
};

inline Value parse_file(const std::string& filepath) {
    std::ifstream f(filepath);
    if (!f.is_open()) throw std::runtime_error("Could not open file: " + filepath);
    std::stringstream buffer;
    buffer << f.rdbuf();
    Parser p;
    return p.parse(buffer.str());
}

} // namespace json
} // namespace fvw

#endif // JSON_UTILS_H
