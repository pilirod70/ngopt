/**
 * @file Log.h
 * @brief Log Functions
 * @author Yu Peng
 * @version 0.18
 * @date 2010-10-26
 */

#ifndef __LOG_H_

#define __LOG_H_

#include "globals.h"

#include <cstdio>

void SetLogFile(std::FILE *newLogFile);
void LogMessage(const char *fmt, ...);
void LogDebug(const char *fmt, ...);
void LogError(const char *fmt, ...);

#endif

