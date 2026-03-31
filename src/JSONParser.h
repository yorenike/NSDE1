#ifndef JSON_PARSER_H
#define JSON_PARSER_H

#include <string>
#include <map>
#include <vector>

class JSONParser {
private:
    std::map<std::string, std::string> data;
    
    // 辅助函数声明
    std::string trim(const std::string& str);
    std::string unquote(const std::string& str);
    std::string removeComments(const std::string& content);
    std::string removeWhitespace(const std::string& str);
    void parse(const std::string& filename);
    void parseObject(const std::string& json, size_t& pos, const std::string& prefix);
    
public:
    // 构造函数
    JSONParser(const std::string& filename);
    
    // 获取值的方法
    std::string getString(const std::string& key) const;
    int getInt(const std::string& key) const;
    double getDouble(const std::string& key) const;
    bool getBool(const std::string& key) const;
    
    // 检查键是否存在
    bool hasKey(const std::string& key) const;
    
    // 带默认值的获取方法
    std::string getStringOrDefault(const std::string& key, const std::string& default_value) const;
    int getIntOrDefault(const std::string& key, int default_value) const;
    double getDoubleOrDefault(const std::string& key, double default_value) const;
    bool getBoolOrDefault(const std::string& key, bool default_value) const;
    
    // 调试工具
    void printAll() const;
    void validateRequiredKeys(const std::vector<std::string>& required_keys) const;
};

#endif // JSON_PARSER_H