/**
 * @file Log.cpp
 * @brief 
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#include "globals.h"
#include "Log.h"

#include <cstdio>
#include <cstdarg>
#include <algorithm>

using namespace std;

static char message[MaxLine];
static FILE *logFile = stderr;

void SetLogFile(FILE *newLogFile)
{
    logFile = newLogFile;
}

void LogMessage(const char *fmt, ...)
{
    if (logFile != NULL)
    {
        va_list ap;
        va_start(ap, fmt);
        vsprintf(message, fmt, ap);
        va_end(ap);

        fprintf(logFile, "Message: %s", message);
    }
}

void LogDebug(const char *fmt, ...)
{
    if (logFile != NULL)
    {
        va_list ap;
        va_start(ap, fmt);
        vsprintf(message, fmt, ap);
        va_end(ap);

        fprintf(logFile, "Debug: %s", message);
    }
}

void LogError(const char *fmt, ...)
{
    if (logFile != NULL)
    {
        va_list ap;
        va_start(ap, fmt);
        vsprintf(message, fmt, ap);
        va_end(ap);

        fprintf(logFile, "Error: %s", message);
    }
}

