
#include "files.h"
#include <algorithm>
#include <filesystem>
#include "common.h"

const char* last_slash(const char* path) {
    const char* result = path;
    while(*path) {
        if(*path == '/') {
            result = path + 1;
        }
        path++;
    }

    return result;
}

bool number_compare(const std::string& a, const std::string& b) {
    int a_number = 0;
    int b_number = 0;

    sscanf(last_slash(a.c_str()), "%d", &a_number);
    sscanf(last_slash(b.c_str()), "%d", &b_number);

    return a_number < b_number;
}

std::vector<std::string> list_files(const char* path, const char* ext) {
    std::vector<std::string> result;
    for(auto iter : std::filesystem::directory_iterator(path)) {
        if(iter.is_directory()) {
            continue;
        }

        if(iter.path().extension() != ext) {
            continue;
        }

        result.push_back(iter.path().string());
    }

    std::sort(result.begin(), result.end(), &number_compare);
    return result;
}

void progress_bar::print_progress(size_t value) {
    int percent = (int)((value * 100) / total);
    if(true) { 
        do_print(percent / 2, 50 - percent / 2, value);
        previous_percent = percent;
    }
}

void progress_bar::done() {
    do_print(50, 0, total);
    printf("\n");
}

progress_bar::progress_bar(size_t total) : total(total) {
    do_print(0, 50, 0);
}

static void format_time(int seconds, char* buffer) {
    int minutes = seconds / 60;
    int hours = minutes / 60;
    int days = hours / 24;

    seconds %= 60;
    minutes %= 60;
    hours %= 24;

    if(days > 0) {
        sprintf(buffer, "%dd %dh %dm %ds", days, hours, minutes, seconds);
    } else if(hours > 0) {
        sprintf(buffer, "%dh %dm %ds", hours, minutes, seconds);
    } else if(minutes > 0) {
        sprintf(buffer, "%dm %ds", minutes, seconds);
    } else if(seconds < 5) {
        sprintf(buffer, "...");
    } else {
        sprintf(buffer, "%ds", seconds);
    }

}
void progress_bar::do_print(size_t markers, size_t spaces, size_t value) {
    printf("\r[");
    for(int i = 0; i < markers; i++) {
        putchar('=');
    }

    for(int i = 0; i < spaces; i++) {
        putchar(' ');
    }

    char buffer[256];
    auto now = std::chrono::system_clock::now();
    auto elapsed = std::chrono::duration_cast<std::chrono::nanoseconds>(now - start).count();

    size_t remaining = 0;

    if(value != 0) {
        remaining = (elapsed * total / value) - elapsed;
    }

    format_time(remaining / 1000000000, buffer);

    printf("] %d%% (%d/%d), %s left          ", (int)((value * 100) / total), (int)value, (int)total, buffer);
    fflush(stdout);
}
