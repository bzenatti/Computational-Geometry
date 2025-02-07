#include <opencv2/core.hpp>
#include <opencv2/imgcodecs.hpp>
#include <opencv2/highgui.hpp>
#include <opencv2/imgproc.hpp>  // For color conversion
#include <iostream>
#include <filesystem>
#include <sys/stat.h>

namespace fs = std::filesystem;
using namespace cv;

bool convert_to_pgm(const std::string& image_path, const std::string& output_path);


int main() {
    std::string image_path = "./baboon.jpg";
    std::string output_path = "./baboon.pgm";

    if (convert_to_pgm(image_path, output_path)) {
        std::cout << "Image conversion successful." << std::endl;
    } else {
        std::cout << "Image conversion failed." << std::endl;
    }

    return 0;
}


bool convert_to_pgm(const std::string& image_path, const std::string& output_path) {
    // Print current working directory
    std::cout << "Current working directory: " << fs::current_path() << std::endl;

    // Check if file exists
    if (!fs::exists(image_path)) {
        std::cout << "File does not exist: " << image_path << std::endl;
        return false;
    } else {
        std::cout << "File exists: " << image_path << std::endl;

        // Check file permissions
        struct stat buffer;
        if (stat(image_path.c_str(), &buffer) == 0) {
            std::cout << "File permissions: " << std::oct << (buffer.st_mode & 0777) << std::endl;
        }
    }

    // Read image in color
    Mat img = imread(image_path, IMREAD_COLOR);
    if (img.empty()) {
        std::cout << "Could not read the image: " << image_path << std::endl;
        return false;
    }

    // Convert to grayscale (PGM format requires grayscale)
    Mat gray_img;
    cvtColor(img, gray_img, COLOR_BGR2GRAY);

    // Save the grayscale image in PGM format
    if (imwrite(output_path, gray_img)) {
        std::cout << "PGM file saved: " << output_path << std::endl;
        return true;
    } else {
        std::cout << "Failed to save PGM file: " << output_path << std::endl;
        return false;
    }
}