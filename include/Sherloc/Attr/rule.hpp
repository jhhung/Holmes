#pragma once

#include <spdlog/fmt/fmt.h>

namespace Sherloc::Attr {

struct Rule{
private:
    int _rule = 0;

    bool _is_enabled = false;
public:
    Rule(int rule, bool is_enable = false): _rule(rule), _is_enabled(is_enable){}
    Rule(const Rule&) = default;
    Rule(Rule&&) = default;

    operator auto() const {
        return _rule;
    }

    void enable(){
        _is_enabled = true;
    }

    void disable(){
        _is_enabled = false;
    }

    [[nodiscard]] inline auto is_enable() const {
        return _is_enabled;
    }

    auto operator<=>(const Rule& rhs) const {
        return _rule <=> rhs._rule;
    }

    [[nodiscard]] auto str() const {
        return fmt::format("{}:{}", _rule, is_enable() ? 'E':'D');
    }
};

}
