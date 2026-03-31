#include "JSONParser.h"
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <cctype>
#include <stdexcept>

// 去除字符串两端的空格
std::string JSONParser::trim(const std::string& str) {
    size_t first = str.find_first_not_of(" \t\n\r");
    if (first == std::string::npos) return "";
    size_t last = str.find_last_not_of(" \t\n\r");
    return str.substr(first, last - first + 1);
}

// 去除字符串两端的引号
std::string JSONParser::unquote(const std::string& str) {
    std::string result = trim(str);
    if (result.size() >= 2 && result.front() == '"' && result.back() == '"') {
        return result.substr(1, result.size() - 2);
    }
    return result;
}

// 去除所有注释（// 开头的行）
std::string JSONParser::removeComments(const std::string& content) {
    std::istringstream iss(content);
    std::string line;
    std::string result;
    
    while (std::getline(iss, line)) {
        size_t comment_pos = line.find("//");
        if (comment_pos != std::string::npos) {
            line = line.substr(0, comment_pos);
        }
        result += line + "\n";
    }
    
    return result;
}

// 去除所有空白字符
std::string JSONParser::removeWhitespace(const std::string& str) {
    std::string result;
    for (char c : str) {
        if (!isspace(c)) {
            result += c;
        }
    }
    return result;
}

// 递归解析对象
void JSONParser::parseObject(const std::string& json, size_t& pos, const std::string& prefix) {
    while (pos < json.size()) {
        // 跳过逗号
        if (json[pos] == ',') {
            pos++;
            continue;
        }
        
        // 遇到结束大括号，返回
        if (json[pos] == '}') {
            pos++;
            return;
        }
        
        // 解析 key
        if (json[pos] != '"') {
            throw std::runtime_error("Expected '\"' at position " + std::to_string(pos));
        }
        pos++; // 跳过开头的 '"'
        
        std::string key;
        while (pos < json.size() && json[pos] != '"') {
            key += json[pos];
            pos++;
        }
        if (json[pos] != '"') {
            throw std::runtime_error("Expected '\"' at position " + std::to_string(pos));
        }
        pos++; // 跳过结尾的 '"'
        
        // 跳过 ':'
        if (json[pos] != ':') {
            throw std::runtime_error("Expected ':' at position " + std::to_string(pos));
        }
        pos++;
        
        // 构建完整路径
        std::string full_key = prefix.empty() ? key : prefix + "." + key;
        
        // 检查 value 的类型
        if (json[pos] == '{') {
            // 对象：递归解析
            pos++; // 跳过 '{'
            parseObject(json, pos, full_key);
        } else if (json[pos] == '"') {
            // 字符串
            pos++; // 跳过开头的 '"'
            std::string value;
            while (pos < json.size() && json[pos] != '"') {
                value += json[pos];
                pos++;
            }
            if (json[pos] != '"') {
                throw std::runtime_error("Expected '\"' at position " + std::to_string(pos));
            }
            pos++; // 跳过结尾的 '"'
            
            data[full_key] = value;
        } else {
            // 数字、布尔值、null
            std::string value;
            while (pos < json.size() && json[pos] != ',' && json[pos] != '}') {
                value += json[pos];
                pos++;
            }
            
            data[full_key] = value;
        }
    }
}

// 主解析函数
void JSONParser::parse(const std::string& filename) {
    // 读取整个文件
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Cannot open file: " + filename);
    }
    
    std::stringstream buffer;
    buffer << file.rdbuf();
    std::string content = buffer.str();
    file.close();
    
    // 去除注释
    content = removeComments(content);
    
    // 去除所有空白字符
    std::string json = removeWhitespace(content);
    
    // 解析 JSON
    size_t pos = 0;
    
    // 跳过开头的 '{'
    if (json[pos] != '{') {
        throw std::runtime_error("JSON must start with '{'");
    }
    pos++;
    
    // 解析顶层对象
    parseObject(json, pos, "");
}

// 构造函数
JSONParser::JSONParser(const std::string& filename) {
    parse(filename);
}

// 获取字符串值
std::string JSONParser::getString(const std::string& key) const {
    auto it = data.find(key);
    if (it == data.end()) {
        throw std::runtime_error("Key not found: " + key);
    }
    return it->second;
}

// 获取整数值
int JSONParser::getInt(const std::string& key) const {
    return std::stoi(getString(key));
}

// 获取浮点数值
double JSONParser::getDouble(const std::string& key) const {
    return std::stod(getString(key));
}

// 获取布尔值
bool JSONParser::getBool(const std::string& key) const {
    std::string val = getString(key);
    std::transform(val.begin(), val.end(), val.begin(), ::tolower);
    return val == "true" || val == "1";
}

// 检查键是否存在
bool JSONParser::hasKey(const std::string& key) const {
    return data.find(key) != data.end();
}

// 带默认值的字符串
std::string JSONParser::getStringOrDefault(const std::string& key, const std::string& default_value) const {
    if (hasKey(key)) {
        return getString(key);
    }
    return default_value;
}

// 带默认值的整数
int JSONParser::getIntOrDefault(const std::string& key, int default_value) const {
    if (hasKey(key)) {
        return getInt(key);
    }
    return default_value;
}

// 带默认值的浮点数
double JSONParser::getDoubleOrDefault(const std::string& key, double default_value) const {
    if (hasKey(key)) {
        return getDouble(key);
    }
    return default_value;
}

// 带默认值的布尔值
bool JSONParser::getBoolOrDefault(const std::string& key, bool default_value) const {
    if (hasKey(key)) {
        return getBool(key);
    }
    return default_value;
}

// 打印所有解析的数据
void JSONParser::printAll() const {
    std::cout << "=== Parsed JSON Data ===" << std::endl;
    for (const auto& pair : data) {
        std::cout << "  " << pair.first << " = " << pair.second << std::endl;
    }
    std::cout << "=========================" << std::endl;
}

// 验证必需的键是否存在
void JSONParser::validateRequiredKeys(const std::vector<std::string>& required_keys) const {
    for (const auto& key : required_keys) {
        if (!hasKey(key)) {
            throw std::runtime_error("Required key missing: " + key);
        }
    }
}