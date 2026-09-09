#pragma once

#include <algorithm>
#include <cctype>
#include <iostream>
#include <string>

#include <sys/wait.h>
#include <cstring>
#include <signal.h>
#include <csignal>
#include <unistd.h>

void trace();

void handler(int signo, siginfo_t* info, void* context);

void warmup_cpptrace();

void segfault_handler_cpptrace ();
