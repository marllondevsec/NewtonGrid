#include <iostream>

int main(int argc, char** argv) {
    std::cout << "gravfield app - use a JSON config file as argument\n";
    if (argc < 2) {
        std::cout << "Usage: " << argv[0] << " examples/configs/single_mass.json\n";
        return 1;
    }
    // TODO: load config, run experiment runner
    return 0;
}
