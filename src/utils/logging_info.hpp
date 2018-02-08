#ifndef LOGGINGINFO_H_
#define LOGGINGINFO_H_
#include <string>
#include <vector>
#include <iostream>
//#include <sstream>

/**@file
 * @brief Logging class for debug
 */

template<typename T>
class LoggingInfo
{
    template<typename A>
    friend std::ostream& operator<<(std::ostream &os, LoggingInfo<A> const& log);
protected:
    std::string     m_name;

    char            m_delimiter;

    int             m_num;

    std::vector<T>  m_buffer;

public:
    LoggingInfo(std::string name, char del = ' ', int num = -1) : 
        m_name(name), m_delimiter(del), m_num(num)
    {}

    virtual ~LoggingInfo() {
        //std::cerr << *this;
    }

    LoggingInfo& operator<<(T const& data) {
        m_buffer.push_back(data);
        return *this;
    }

    int length() {
        return m_buffer.size();
    }
};

template<typename T>
std::ostream& operator<<(std::ostream &os, LoggingInfo<T> const& log)
{
    for (int i = 0; i < (int)log.m_buffer.size(); i++) {
        if (log.m_num < 0) {
            if (log.m_delimiter == '\n')
                os << "# " << i + 1 << ": " << log.m_buffer[i] << log.m_delimiter;
            else 
                os << log.m_buffer[i] << log.m_delimiter;
        } else {
            if (i % log.m_num == 0)
                os << "# " << i + 1 << " - " << i + log.m_num << ": " 
                    << log.m_buffer[i] << log.m_delimiter;
            else
                os << log.m_buffer[i] << log.m_delimiter;
        }
        if (log.m_num > 0 && (i + 1) % log.m_num == 0)
            os << std::endl;
    }
    if (log.m_buffer.size() != 0) os << std::endl;

    return os;
}

#endif
