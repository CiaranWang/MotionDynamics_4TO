#include "md.h"
#include "ReadInput.h"

#include <unordered_map>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <vector>
#include <fstream>
#include <string>

constexpr double sigma = 100;

// ----------------------
// Gaussian blur for density map
// ----------------------
std::vector<std::vector<double>> gaussianBlur2D(
    const std::vector<std::vector<double>>& input, double sigma)
{
    int nx = input.size();
    int ny = input[0].size();
    std::vector<std::vector<double>> temp(nx, std::vector<double>(ny, 0.0));
    std::vector<std::vector<double>> output(nx, std::vector<double>(ny, 0.0));

    int radius = static_cast<int>(3.0 * sigma);
    std::vector<double> kernel(2 * radius + 1);
    double sum = 0.0;
    for (int i = -radius; i <= radius; i++) {
        kernel[i + radius] = std::exp(-(i * i) / (2 * sigma * sigma));
        sum += kernel[i + radius];
    }
    for (auto& k : kernel) k /= sum;

    // horizontal blur
    for (int x = 0; x < nx; x++)
        for (int y = 0; y < ny; y++) {
            double val = 0.0;
            for (int k = -radius; k <= radius; k++) {
                int ix = std::min(std::max(x + k, 0), nx - 1);
                val += input[ix][y] * kernel[k + radius];
            }
            temp[x][y] = val;
        }

    // vertical blur
    for (int x = 0; x < nx; x++)
        for (int y = 0; y < ny; y++) {
            double val = 0.0;
            for (int k = -radius; k <= radius; k++) {
                int iy = std::min(std::max(y + k, 0), ny - 1);
                val += temp[x][iy] * kernel[k + radius];
            }
            output[x][y] = val;
        }

    return output;
}

// ----------------------
// Create density map
// ----------------------
std::vector<std::vector<double>> createDensityMap(
    const std::vector<Detection>& detections,
    int grid_x, int grid_y,
    double pen_width, double pen_height)
{
    std::vector<std::vector<double>> density(grid_x, std::vector<double>(grid_y, 0.0));
    double dx = pen_width / grid_x;
    double dy = pen_height / grid_y;

    for (auto& det : detections) {
        int ix = static_cast<int>(det.cen_x / dx);
        int iy = static_cast<int>(det.cen_y / dy);
        if (ix >= 0 && ix < grid_x && iy >= 0 && iy < grid_y)
            density[ix][iy] += 1.0;
    }

    return gaussianBlur2D(density, sigma);
}

// ----------------------
// Accumulated density ahead
// ----------------------
double accumulatedDensityAhead(
    const Detection& det,
    const std::vector<std::vector<double>>& density,
    double dx, double dy)
{
    int nx = density.size();
    int ny = density[0].size();
    double acc = 0.0;
    double factor_sum = 0.0;

    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            double x_rel = (i + 0.5) * dx - det.cen_x;
            double y_rel = (j + 0.5) * dy - det.cen_y;
            double r = std::sqrt(x_rel * x_rel + y_rel * y_rel);
            if (r > 300.0) continue;

            double theta = std::atan2(y_rel, x_rel) - det.angle;
            while (theta < -kPI) theta += 2 * kPI;
            while (theta > kPI) theta -= 2 * kPI;
            if (theta<-kPI / 4 || theta>kPI / 4) continue;

            double weight = 1.0 - r / 150.0;
            acc += weight * density[i][j];
            factor_sum += weight;
        }
    }
    if (factor_sum > 0) acc /= factor_sum;
    return acc;
}

// ----------------------
// Build segments with frame_window and interpolate missing frames
// ----------------------
std::vector<Segment> buildSegments(
    std::vector<Detection>& detections, long frame_window,
    const std::vector<std::vector<double>>& density_map,
    double dx, double dy)
{
    std::vector<Segment> segments;
    if (detections.empty()) return segments;

    std::sort(detections.begin(), detections.end(), [](const Detection& a, const Detection& b) { return a.custom_frame < b.custom_frame; });

    Segment seg;
    Detection* prev = nullptr;

    for (size_t i = 0; i < detections.size(); i++) {
        Detection& det = detections[i];
        // compute local density and ahead_density
        det.local_density = density_map[std::min(std::max(int(det.cen_x / dx), 0), int(density_map.size() - 1))]
            [std::min(std::max(int(det.cen_y / dy), 0), int(density_map[0].size() - 1))];
        det.ahead_density = accumulatedDensityAhead(det, density_map, dx, dy);

        if (!prev) {
            seg.points.push_back(det);
            prev = &seg.points.back();
            continue;
        }

        long gap = det.custom_frame - prev->custom_frame;
        if (gap > frame_window) {
            // save previous segment
            segments.push_back(seg);
            seg.points.clear();
            seg.points.push_back(det);
            prev = &seg.points.back();
            continue;
        }

        // interpolate missing frames
        for (long f = 1; f < gap; f++) {
            double alpha = double(f) / gap;
            Detection interp;
            interp.ID = det.ID;
            interp.custom_frame = prev->custom_frame + f;

            interp.cen_x = prev->cen_x + alpha * (det.cen_x - prev->cen_x);
            interp.cen_y = prev->cen_y + alpha * (det.cen_y - prev->cen_y);
            interp.dir_x = prev->dir_x + alpha * (det.dir_x - prev->dir_x);
            interp.dir_y = prev->dir_y + alpha * (det.dir_y - prev->dir_y);
            interp.angle = std::atan2(interp.dir_y, interp.dir_x);
            if (interp.angle < 0) interp.angle += 2.0 * kPI;
            interp.interpolated = true;

            // densities
            interp.local_density = density_map[std::min(std::max(int(interp.cen_x / dx), 0), int(density_map.size() - 1))]
                [std::min(std::max(int(interp.cen_y / dy), 0), int(density_map[0].size() - 1))];
            interp.ahead_density = accumulatedDensityAhead(interp, density_map, dx, dy);

            seg.points.push_back(interp);
        }

        // distance, speed, angle change
        det.distance_from_last = std::sqrt((det.cen_x - prev->cen_x) * (det.cen_x - prev->cen_x) +
            (det.cen_y - prev->cen_y) * (det.cen_y - prev->cen_y));
        det.speed = det.distance_from_last;
        det.acceleration = det.speed - prev->speed;
        double dtheta = det.angle - prev->angle;
        while (dtheta < -kPI) dtheta += 2 * kPI;
        while (dtheta > kPI) dtheta -= 2 * kPI;
        det.angle_change = dtheta;

        seg.points.push_back(det);
        prev = &seg.points.back();
    }

    if (!seg.points.empty()) segments.push_back(seg);
    return segments;
}

// ----------------------
// Simple moving average smoothing
// ----------------------
void smoothVector(std::vector<double>& data, int window = 5) {
    if (data.size() < 2) return;
    std::vector<double> temp = data;
    int half = window / 2;
    for (size_t i = 0; i < data.size(); i++) {
        double sum = 0.0;
        int count = 0;
        for (int j = std::max<int>(0, i - half); j <= std::min<int>(data.size() - 1, i + half); j++) {
            sum += temp[j];
            count++;
        }
        data[i] = sum / count;
    }
}

// ----------------------
// Main engine
// ----------------------
int motion_dynamics_run(const std::string& input,
    const std::string& output_prefix,
    long frame_window,
    bool smooth)
{
    auto detections_all = ReadInput_Tzayhri(input);
    std::cout << "Read " << detections_all.size() << " detections from " << input << "\n";

    const double pen_width = 2560.0;
    const double pen_height = 1440.0;
    const int grid_x = 256;
    const int grid_y = 144;
    double dx = pen_width / grid_x;
    double dy = pen_height / grid_y;

    auto density_map = createDensityMap(detections_all, grid_x, grid_y, pen_width, pen_height);

    std::unordered_map<int, std::vector<Detection>> per_id;
    for (auto& d : detections_all) per_id[d.ID].push_back(d);

    for (auto& [id, dets] : per_id) {
        auto segments = buildSegments(dets, frame_window, density_map, dx, dy);

        // optional smoothing
        if (smooth) {
            for (auto& seg : segments) {
                std::vector<double> speeds, accs, angles;
                for (auto& det : seg.points) {
                    speeds.push_back(det.speed);
                    accs.push_back(det.acceleration);
                    angles.push_back(det.angle_change);
                }
                smoothVector(speeds, 5);
                smoothVector(accs, 5);
                smoothVector(angles, 5);
                for (size_t k = 0; k < seg.points.size(); k++) {
                    seg.points[k].speed = speeds[k];
                    seg.points[k].acceleration = accs[k];
                    seg.points[k].angle_change = angles[k];
                }
            }
        }

        std::string filename = output_prefix + "_ID" + std::to_string(id) + ".csv";
        std::ofstream fout(filename);
        fout << "segment_index,frame,cen_x,cen_y,dir_x,dir_y,angle,ahead_density,local_density,speed,acceleration,angle_change,distance_from_last,interpolated\n";

        for (size_t s = 0; s < segments.size(); s++) {
            for (auto& det : segments[s].points) {
                fout << s << "," << det.custom_frame << "," << det.cen_x << "," << det.cen_y
                    << "," << det.dir_x << "," << det.dir_y << "," << det.angle
                    << "," << det.ahead_density << "," << det.local_density
                    << "," << det.speed << "," << det.acceleration
                    << "," << det.angle_change << "," << det.distance_from_last
                    << "," << (det.interpolated ? "true" : "false") << "\n";
            }
        }
        std::cout << "Written " << filename << " with " << segments.size() << " segments\n";
    }

    return 0;
}