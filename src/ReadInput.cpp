#include "ReadInput.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <cmath>
#include <vector>

// minimal CSV splitter
static std::vector<std::string> splitCSV(const std::string& line) {
    std::vector<std::string> out;
    std::string cur;
    bool in_quotes = false;

    for (char c : line) {
        if (c == '"') {
            in_quotes = !in_quotes;
        }
        else if (c == ',' && !in_quotes) {
            out.push_back(cur);
            cur.clear();
        }
        else {
            cur.push_back(c);
        }
    }
    out.push_back(cur);
    return out;
}

// Parser for coord_paper4.csv
std::vector<Detection> ReadInput_Tzayhri(const std::string& filename) {
    std::vector<Detection> out;
    std::ifstream fin(filename);

    if (!fin.is_open()) {
        std::cerr << "Error opening: " << filename << "\n";
        return out;
    }

    std::string line;
    std::getline(fin, line);  // skip header

    while (std::getline(fin, line)) {
        auto cols = splitCSV(line);
        if (cols.size() < 31) continue;

        try {
            Detection d;
            d.ID = std::stoi(cols[2]);
            d.custom_frame = std::stol(cols[18]);
            d.cen_x = std::stod(cols[29]);
            d.cen_y = std::stod(cols[30]);

            double tlx = std::stod(cols[3]);
            double tly = std::stod(cols[4]);
            double trx = std::stod(cols[5]);
            double try_ = std::stod(cols[6]);
            double brx = std::stod(cols[7]);
            double bry = std::stod(cols[8]);
            double blx = std::stod(cols[9]);
            double bly = std::stod(cols[10]);

            d.dir_x = (tlx + trx) - (brx + blx);
            d.dir_y = (tly + try_) - (bry + bly);

            d.angle = std::atan2(d.dir_y, d.dir_x);
            if (d.angle < 0) d.angle += 2.0 * kPI;

            out.push_back(d);
        }
        catch (...) { continue; }
    }
    return out;
}