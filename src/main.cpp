#include "md.h"
#include <iostream>
#include <string>
#include <cstdlib>

#ifdef _WIN32
#include <windows.h>
#else
#include <unistd.h>
#include <libgen.h>
#include <limits.h>
#include <getopt.h>
#endif

#ifdef _WIN32
int optind = 1;
char* optarg = nullptr;
int getopt_simple(int argc, char* argv[], const char* opts) {
    if (optind >= argc) return -1;
    char* arg = argv[optind];
    if (arg[0] != '-') return -1;
    char c = arg[1];
    const char* p = strchr(opts, c);
    if (!p) return '?';
    if (*(p + 1) == ':') {
        if (optind + 1 >= argc) return '?';
        optarg = argv[++optind];
    }
    optind++;
    return c;
}
#endif

static std::string getExecutableDir() {
#ifdef _WIN32
    char buffer[MAX_PATH];
    GetModuleFileNameA(NULL, buffer, MAX_PATH);
    std::string fullpath(buffer);
    size_t pos = fullpath.find_last_of("\\/");
    if (pos != std::string::npos)
        return fullpath.substr(0, pos);
    return ".";
#else
    char buff[PATH_MAX];
    ssize_t len = readlink("/proc/self/exe", buff, sizeof(buff) - 1);
    if (len != -1) { buff[len] = '\0'; return std::string(dirname(buff)); }
    return ".";
#endif
}

int main(int argc, char* argv[]) {

    std::string exec_dir = getExecutableDir();
    std::string input_file;
    std::string output_prefix = "results";
    long frame_window = 200;
    bool smooth_flag = false;

    // Handle global flags --version and --update
    for (int i = 1; i < argc; i++) {
        std::string arg = argv[i];
        if (arg == "--version" || arg == "-V") {
            std::cout << "MotionDynamics version " << MD_VERSION << "\n";
            return 0;
        }
        if (arg == "--update" || arg == "-U") {
#ifdef _WIN32
            std::cout << "Auto-update not supported on Windows.\n";
            return 1;
#else
            std::cout << "Updating MotionDynamics...\n";
            std::cout << "Executable directory: " << exec_dir << "\n";
            std::string cmd =
                "cd \"" + exec_dir + "/..\" && "
                "git pull && "
                "cmake -S . -B build && "
                "cmake --build build";
            int ret = std::system(cmd.c_str());
            if (ret == 0) std::cout << "Update & rebuild complete.\n";
            else std::cout << "Update failed.\n";
            return ret;        //heya
#endif
        }
        if (arg == "--smooth") smooth_flag = true; // optional smoothing
    }

#ifndef _WIN32
    static struct option long_opts[] = { {"frame_window", required_argument, 0, 'f'}, {0,0,0,0} };
    int opt;
    while ((opt = getopt_long(argc, argv, "i:o:f:", long_opts, NULL)) != -1) {
#else
    int opt;
    while ((opt = getopt_simple(argc, argv, "i:o:f:")) != -1) {
#endif
        switch (opt) {
        case 'i': input_file = optarg; break;
        case 'o': output_prefix = optarg; break;
        case 'f': frame_window = std::stol(optarg); break;
        default:
            std::cerr << "Usage:\n"
                << "  MotionDynamics -i file.csv -o out --frame_window N [--smooth]\n"
                << "  MotionDynamics --update | -U\n"
                << "  MotionDynamics --version | -V\n";
            return 1;
        }
    }

    if (input_file.empty()) {
        std::cerr << "Error: no input file (-i)\n";
        return 1;
    }

    return motion_dynamics_run(input_file, output_prefix, frame_window, smooth_flag);
    }