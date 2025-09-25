#include "logger.h"

Logger* Logger::instance = nullptr; 
std::string Logger::logfilepath = "log.txt";
std::ofstream Logger::logfile;
uint Logger::layer = 0;
VERBOSITY Logger::verbosity = VERBOSITY::TWO;

Logger::Logger() {}
Logger::~Logger() {
    logfile.close();
}

void Logger::log(std::string message, LOGTYPE type, VERBOSITY verbosity_) {
    if(verbosity_ > Logger::verbosity && type != LOGTYPE::WARNING) return;
    if(not logfile.is_open()) {
        logfile.open(Logger::logfilepath,std::ios::app);
    }
    std::time_t t = std::time(nullptr);
    std::tm tm = *std::localtime(&t);
    std::stringstream timestamp_sstr;
    timestamp_sstr << std::put_time(&tm, "[%Y-%m-%d_%H-%M-%S]");
    std::string timestamp = timestamp_sstr.str();
    
    std::string typeIndicator;
    switch (type)
    {
    case LOGTYPE::WARNING:
        typeIndicator = "[WARNING]";
        break;
    case LOGTYPE::INFO:
        typeIndicator = "[INFO]";
        break;
    case LOGTYPE::DEBUG:
        typeIndicator = "[DEBUG]";
        break;    
    default:
        break;
    }

    std::string additionalTabs = "\t\t";
    for(uint i = 0; i < Logger::layer; i++) additionalTabs += "\t";

    logfile << timestamp << ":\t" << typeIndicator << additionalTabs << message << std::endl;
}

void Logger::setLogfilepath(std::string logfilepath) {
    logfile.close();
    this->logfilepath = logfilepath;
}

void Logger::clearLogfile() {
    logfile.close();
    logfile.open(Logger::logfilepath,std::ios::trunc);
}

Logger& Logger::getInstance() {
    // If the instance doesn't exist, create it
    if (!instance) {
        instance = new Logger();
    }
    return *instance;
}
