#include "logger.h"
#include <iostream>
#include <ctime>
#include <cstdarg>
#include <cstring>

using namespace std;

Logger::Logger(int log_level) : log_level(static_cast<LogLevel>(log_level)) {}

void Logger::debugl2(const char *format, ...) {
    va_list args;
    va_start(args, format);
    log(formatString(format, args), DEBUGL2);
    va_end(args);
}

void Logger::debug(const char *const format, ...) {
    va_list args;
    va_start(args, format);
    log(formatString(format, args), DEBUG);
    va_end(args);
}

void Logger::info(const char *const format, ...) {
    va_list args;
    va_start(args, format);
    log(formatString(format, args), INFO);
    va_end(args);
}

void Logger::warn(const char *const format, ...) {
    va_list args;
    va_start(args, format);
    log(formatString(format, args), WARN);
    va_end(args);

}

void Logger::error(const char *const format, ...) {
    va_list args;
    va_start(args, format);
    log(formatString(format, args), ERROR);
    va_end(args);
}

void Logger::fatal(const char *const format, ...) {
    va_list args;
    va_start(args, format);
    log(formatString(format, args), FATAL);
    va_end(args);
}

void Logger::log(const string &s, LogLevel log_level) {
    if (log_level > this->log_level)
        return;
    time_t ctt = time(nullptr);
    char *time = asctime(localtime(&ctt));
    time[strlen(time) - 1] = '\0';
    const char *s_c_str = s.c_str();
    mtx.lock();
    printf("%s: %s\n", time, s_c_str);
//    cout << time << ": " << s << endl;
    mtx.unlock();
}

string Logger::formatString(const char *const format, va_list args) {
    char buffer[FORMAT_LENGTH];
    vsnprintf(buffer, FORMAT_LENGTH, format, args);
    return string(buffer);
}

string Logger::formatString(const char *const format, ...) {
    va_list args;
    va_start(args, format);
    string s = formatString(format, args);
    va_end(args);
    return s;
}

