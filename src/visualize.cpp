#include "visualize.hpp"
#include <iostream>
#include <unordered_set>
#include <cstdio>
#include <unordered_map>
#include <cmath>

// 在控制台打印网格、起点、终点和路径
void print_grid_with_path(
    const std::vector<std::vector<int>>& grid,
    Point start, Point goal,
    const std::vector<Point>& path
){
    int H = grid.size();
    int W = grid[0].size();
    std::vector<std::string> out(H, std::string(W, '?'));

    for (int y = 0; y < H; ++y)
        for (int x = 0; x < W; ++x)
            out[y][x] = (grid[y][x] == 1 ? '#' : '.');

    for (auto& p : path)
        out[p.y][p.x] = '*';

    out[start.y][start.x] = 'S';
    out[goal.y][goal.x]   = 'G';

    for (int y = 0; y < H; ++y)
        std::cout << out[y] << "\n";
}

// 将网格、起点、终点和路径写入 PPM 文件，支持缩放
void write_ppm_scaled(
    const std::string &filename,
    const std::vector<std::vector<int>>& grid,
    Point start, Point goal,
    const std::vector<Point>& path,
    int scale
){
    int H = grid.size();
    int W = grid[0].size();

    int outW = W * scale;
    int outH = H * scale;

    FILE *f = fopen(filename.c_str(), "wb");
    if (!f) {
        std::cerr << "无法创建 PPM 文件\n";
        return;
    }

    fprintf(f, "P6\n%d %d\n255\n", outW, outH);

    std::unordered_set<long long> pathset;
    for (auto &p : path)
        pathset.insert(((long long)p.x << 32) | p.y);

    long long startkey = ((long long)start.x << 32) | start.y;
    long long goalkey  = ((long long)goal.x << 32) | goal.y;

    unsigned char pixel[3];

    for (int y = 0; y < H; ++y) {
        for (int sy = 0; sy < scale; ++sy) {
            for (int x = 0; x < W; ++x) {

                long long key = ((long long)x << 32) | y;

                if (key == startkey) {
                    pixel[0]=0; pixel[1]=255; pixel[2]=0;
                } else if (key == goalkey) {
                    pixel[0]=0; pixel[1]=0; pixel[2]=255;
                } else if (pathset.count(key)) {
                    pixel[0]=255; pixel[1]=0; pixel[2]=0;
                } else if (grid[y][x] == 1) {
                    pixel[0]=0; pixel[1]=0; pixel[2]=0;
                } else {
                    pixel[0]=255; pixel[1]=255; pixel[2]=255;
                }

                for (int sx = 0; sx < scale; ++sx)
                    fwrite(pixel, 1, 3, f);
            }
        }
    }
    fclose(f);
}

void write_ppm_segments(
    const std::string& filename,
    const std::vector<std::vector<int>>& grid,
    const std::vector<Point>& path,
    const std::vector<Segment>& segments,
    int scale
){
    int H = grid.size();
    int W = grid[0].size();
    int outW = W * scale;
    int outH = H * scale;

    FILE* f = fopen(filename.c_str(), "wb");
    if (!f) {
        std::cerr << "无法创建分段可视化文件\n";
        return;
    }

    fprintf(f, "P6\n%d %d\n255\n", outW, outH);

    // —— 构造路径点 → segment编号 的映射 ——
    std::unordered_map<long long, int> seg_map;
    for (int sid = 0; sid < segments.size(); ++sid) {
        Point s = segments[sid].start;
        Point e = segments[sid].end;

        bool inside = false;
        for (const auto& p : path) {
            if (p.x == s.x && p.y == s.y) inside = true;
            if (inside) {
                long long key = ((long long)p.x << 32) | p.y;
                seg_map[key] = sid;

                if (p.x == e.x && p.y == e.y) break;
            }
        }
    }

    // 一套循环颜色
    const unsigned char colors[][3] = {
        {255,0,0}, {0,255,0}, {0,0,255},
        {255,255,0}, {255,0,255}, {0,255,255},
        {255,128,0}, {128,0,255}
    };
    int color_count = sizeof(colors)/sizeof(colors[0]);

    unsigned char pixel[3];

    for (int y = 0; y < H; ++y) {
        for (int sy = 0; sy < scale; ++sy) {

            for (int x = 0; x < W; ++x) {
                long long key = ((long long)x << 32) | y;

                if (seg_map.count(key)) {
                    int sid = seg_map[key] % color_count;
                    pixel[0] = colors[sid][0];
                    pixel[1] = colors[sid][1];
                    pixel[2] = colors[sid][2];
                }
                else if (grid[y][x] == 1) {
                    pixel[0] = 0; pixel[1] = 0; pixel[2] = 0;
                } else {
                    pixel[0] = 255; pixel[1] = 255; pixel[2] = 255;
                }

                for (int sx = 0; sx < scale; ++sx)
                    fwrite(pixel, 1, 3, f);
            }
        }
    }

    fclose(f);
}

void draw_polygon_ppm(
    const std::string& file,
    const Polygon& poly,
    const std::vector<std::vector<int>>& grid,
    int scale)
{
    int H = grid.size(); 
    int W = grid[0].size();

    int outW = W * scale;
    int outH = H * scale;

    FILE* f = fopen(file.c_str(), "wb");
    fprintf(f, "P6\n%d %d\n255\n", outW, outH);

    unsigned char px[3];

    bool isRect = (poly.size() == 4);

    for (int y = 0; y < H; ++y) {
        for (int sy = 0; sy < scale; ++sy) {

            for (int x = 0; x < W; ++x) {

                bool border = false;

                if (isRect) {
                    // -----------------------------
                    // 🎯 矩形 corridor 精准判定边界
                    // -----------------------------
                    for (int i = 0; i < 4; ++i) {
                        P2 a = poly[i];
                        P2 b = poly[(i + 1) % 4];

                        // 边向量
                        double vx = b.x - a.x;
                        double vy = b.y - a.y;

                        // cell 中心坐标
                        double cx = x + 0.5;
                        double cy = y + 0.5;

                        // 点到边向量
                        double wx = cx - a.x;
                        double wy = cy - a.y;

                        double L = std::sqrt(vx * vx + vy * vy) + 1e-6;
                        double dist = std::fabs(vx * wy - vy * wx) / L;

                        // 判断投影是否在线段范围内
                        double t = (vx * wx + vy * wy) / (L * L);
                        if (t >= 0 && t <= 1 && dist < 0.2) {
                            border = true;
                            break;
                        }
                    }
                }
                else {
                    // -----------------------------
                    // 原有多边形边界检测逻辑
                    // -----------------------------
                    for(int i=0; i<poly.size(); ++i){
                        P2 a=poly[i];
                        P2 b=poly[(i+1)%poly.size()];

                        double vx = b.x - a.x;
                        double vy = b.y - a.y;

                        double cx = x + 0.5;
                        double cy = y + 0.5;

                        double wx = cx - a.x;
                        double wy = cy - a.y;

                        double cross = std::fabs(vx*wy - vy*wx) /
                                    (std::sqrt(vx*vx + vy*vy) + 1e-6);

                        if (cross < 0.1) {
                            border = true;
                            break;
                        }
                    }
                }

                // -----------------------------
                // 上色策略与原版本完全一致
                // -----------------------------
                if(grid[y][x]==1){
                    px[0]=0; px[1]=0; px[2]=0;
                }
                else if(border){
                    px[0]=255; px[1]=255; px[2]=0;  // yellow corridor boundary
                }
                else{
                    px[0]=255; px[1]=255; px[2]=255;
                }

                for (int sx = 0; sx < scale; ++sx)
                    fwrite(px,1,3,f);
            }
        }
    }
    fclose(f);
}

