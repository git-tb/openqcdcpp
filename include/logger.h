#ifndef LOGGER_H
#define LOGGER_H

#include <string>
#include <fstream>
#include <ctime>   // std::time, std::local_time
#include <iomanip> // std::put_time

enum class LOGTYPE { WARNING, INFO, DEBUG};
enum class VERBOSITY { ZERO, ONE, TWO };

class Logger {
    public:
        // Static method to access the singleton instance
        static Logger& getInstance();
        void setLogfilepath(std::string logfilepath);
        void clearLogfile();
        void log(std::string message, LOGTYPE type = LOGTYPE::INFO, VERBOSITY verbosity = VERBOSITY::ONE);
    
        // Delete the copy constructor and assignment operator
        Logger(const Logger&) = delete;
        Logger& operator=(const Logger&) = delete;
        
        static unsigned int layer;
        static VERBOSITY verbosity;
        
    private:
        Logger();
        ~Logger();
    
        static Logger* instance;
        static std::string logfilepath;
        static std::ofstream logfile;
};

#endif