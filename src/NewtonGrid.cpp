// main.cpp
#include <iostream>
#include <exception>
#include <string>
#include <vector>
#include <tuple>
#include <random>
#include <thread>
#include <chrono>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <cmath>

// Inclui os cabeçalhos do NewtonGrid com caminhos relativos
#include "data/Grid.hpp"
#include "data/FieldData.hpp"
#include "physics/GravitationalField.hpp"
#include "io/VTIWriter.hpp"

// ANSI color helpers (funciona em Linux/macOS; tenta habilitar VT em Windows quando disponível)
#define ANSI_RED       "\033[31m"
#define ANSI_RED_BOLD  "\033[1;31m"
#define ANSI_RESET     "\033[0m"

#ifdef _WIN32
  #include <windows.h>
  static void enable_virtual_terminal() {
      HANDLE hOut = GetStdHandle(STD_OUTPUT_HANDLE);
      if (hOut == INVALID_HANDLE_VALUE) return;
      DWORD dwMode = 0;
      if (!GetConsoleMode(hOut, &dwMode)) return;
      dwMode |= ENABLE_VIRTUAL_TERMINAL_PROCESSING;
      SetConsoleMode(hOut, dwMode);
  }
#else
  static void enable_virtual_terminal() { /* noop em terminais POSIX */ }
#endif

using namespace NewtonGrid;

// ----------------------------- util helpers -----------------------------
static void print_banner_animated() {
    // ASCII banner (pode ajustar). Será impresso de cima para baixo com pequeno delay.
    const std::vector<std::string> banner = {
R"(   .+------+     +------+     +------+ )",
R"( .' |    .'|    /|     /|     |      | )",
R"(+---+--+'  |   +-+----+ |     +------+ )",
R"(|   |  |   |   | |    | |     |      | )",
R"(|  ,+--+---+   | +----+-+     +------+ )",
R"(|.'    | .'    |/     |/      |      | )",
R"(+------+'      +------+       +------+ )",
"",
R"( _______                              _______       _     _ )",
R"((_______)              _             (_______)     (_)   | |)",
R"( _     _ _____ _ _ _ _| |_ ___  ____  _   ___  ____ _  __| |)",
R"(| |   | | ___ | | | (_   _) _ \|  _ \| | (_  |/ ___) |/ _  |)",
R"(| |   | | ____| | | | | || |_| | | | | |___) | |   | ( (_| |)",
R"(|_|   |_|_____)\___/   \__)___/|_| |_|\_____/|_|   |_|\____|)",
R"(                                               By MarllonDevSec)"
    };

    // Ativa cor (usar ANSI_RED_BOLD para vermelho em negrito)
    std::cout << ANSI_RED_BOLD;

    for (const auto &line : banner) {
        std::cout << line << std::endl;
        std::this_thread::sleep_for(std::chrono::milliseconds(60));
    }

    // Reseta cor e espaçamento final
    std::cout << ANSI_RESET << std::endl;
    std::this_thread::sleep_for(std::chrono::milliseconds(120));
}

static std::size_t prompt_size_t(const std::string &msg, std::size_t default_val) {
    std::cout << msg << " [" << default_val << "]: ";
    std::string s; std::getline(std::cin, s);
    if (s.empty()) return default_val;
    try { return static_cast<std::size_t>(std::stoul(s)); }
    catch (...) { return default_val; }
}

static double prompt_double(const std::string &msg, double default_val) {
    std::cout << msg << " [" << default_val << "]: ";
    std::string s; std::getline(std::cin, s);
    if (s.empty()) return default_val;
    try { return std::stod(s); }
    catch (...) { return default_val; }
}

static int prompt_int(const std::string &msg, int default_val) {
    std::cout << msg << " [" << default_val << "]: ";
    std::string s; std::getline(std::cin, s);
    if (s.empty()) return default_val;
    try { return std::stoi(s); }
    catch (...) { return default_val; }
}

// Distance between two points
static double dist3(double x1,double y1,double z1,double x2,double y2,double z2){
    const double dx = x1-x2, dy = y1-y2, dz = z1-z2;
    return std::sqrt(dx*dx + dy*dy + dz*dz);
}

// ----------------------------- main -----------------------------
int main() {
    // habilita processamento de sequências ANSI no Windows quando possível
    enable_virtual_terminal();

    print_banner_animated();

    std::cout << "========================================\n";
    std::cout << "   NewtonGrid: Simulação Gravitacional 3D\n";
    std::cout << "========================================\n\n";

    try {
        // ---------------- menu de presets / customização do grid ----------------
        std::cout << "Escolha um preset de grid ou personalize:\n";
        std::cout << "  1) Pequeno (32x32x32, dx=0.1)\n";
        std::cout << "  2) Médio  (50x50x50, dx=0.1)\n";
        std::cout << "  3) Grande (100x100x100, dx=0.1)\n";
        std::cout << "  4) Personalizado\n";
        int grid_choice = prompt_int("Opção", 2);

        std::size_t nx=50, ny=50, nz=50;
        double dx=0.1, dy=0.1, dz=0.1;

        if (grid_choice == 1) {
            nx = ny = nz = 32; dx = dy = dz = 0.1;
        } else if (grid_choice == 2) {
            nx = ny = nz = 50; dx = dy = dz = 0.1;
        } else if (grid_choice == 3) {
            nx = ny = nz = 100; dx = dy = dz = 0.1;
        } else {
            nx = prompt_size_t("nx (pontos em X)", 50);
            ny = prompt_size_t("ny (pontos em Y)", 50);
            nz = prompt_size_t("nz (pontos em Z)", 50);
            dx = prompt_double("dx (espaçamento X)", 0.1);
            dy = prompt_double("dy (espaçamento Y)", 0.1);
            dz = prompt_double("dz (espaçamento Z)", 0.1);
        }

        // Simple safety checks
        if (nx == 0 || ny == 0 || nz == 0) {
            std::cerr << "Dimensões inválidas\n"; return 1;
        }
        if (dx <= 0 || dy <= 0 || dz <= 0) {
            std::cerr << "Espaçamento inválido\n"; return 1;
        }

        const double ox = -(static_cast<double>(nx) * dx) / 2.0;
        const double oy = -(static_cast<double>(ny) * dy) / 2.0;
        const double oz = -(static_cast<double>(nz) * dz) / 2.0;

        // Create Grid
        std::cout << "\n1. Criando grid cartesiano 3D...\n";
        Grid grid({nx, ny, nz}, {dx, dy, dz}, {ox, oy, oz});
        const std::size_t expected_points = nx * ny * nz;
        if (grid.total_points() != expected_points) {
            std::cerr << "ERRO: grid.total_points() mismatch\n"; return 1;
        }

        std::cout << "   • Dimensões: " << nx << " × " << ny << " × " << nz << " pontos\n";
        std::cout << "   • Espaçamento: " << dx << " × " << dy << " × " << dz << "\n";
        std::cout << "   • Origem: (" << ox << ", " << oy << ", " << oz << ")\n";
        std::cout << "   • Pontos totais: " << grid.total_points() << "\n";
        std::cout << "   • Centro físico: (" << grid.center()[0] << ", "
                  << grid.center()[1] << ", " << grid.center()[2] << ")\n";
        std::cout << "   ✓ Grid criado com sucesso!\n\n";

        // ---------------- FieldData allocation ----------------
        std::cout << "2. Alocando campo vetorial...\n";
        FieldData field(grid, true);
        std::cout << "   • FieldData alocado com " << field.size() << " pontos\n\n";

        // ---------------- physics parameters ----------------
        std::cout << "3. Parâmetros físicos (pressione Enter para manter valor padrão)\n";
        double G = prompt_double("Constante gravitacional G", 1.0);
        double softening = prompt_double("Softening (epsilon)", 1e-3);

        if (G <= 0) { std::cerr << "G deve ser > 0\n"; return 1; }
        if (softening < 0) { std::cerr << "Softening deve ser >= 0\n"; return 1; }

        GravitationalField gravity(G, softening);
        std::cout << "   • G = " << gravity.gravitational_constant()
                  << ", softening = " << gravity.softening() << "\n\n";

        // ---------------- massa input menu ----------------
        std::cout << "4. Como deseja definir as massas?\n";
        std::cout << "  1) Manual (digitar posições)\n";
        std::cout << "  2) Simétrico (pares ±x)\n";
        std::cout << "  3) Aleatório com distância mínima\n";
        std::cout << "  4) Ler de arquivo (x y z mass por linha)\n";
        int mass_mode = prompt_int("Opção", 2);

        std::size_t n_masses = prompt_size_t("Número total de massas", 2);

        std::vector<std::tuple<double,double,double,double>> masses; // x,y,z,m

        if (mass_mode == 1) {
            std::cout << "Entrada manual: digite x y z mass (uma por linha)\n";
            for (std::size_t i = 0; i < n_masses; ++i) {
                std::cout << "Massa " << (i+1) << " > ";
                std::string line; std::getline(std::cin, line);
                if (line.empty()) { --i; continue; }
                std::istringstream iss(line);
                double x,y,z,m; if (!(iss >> x >> y >> z >> m)) {
                    std::cout << "Entrada inválida, tente novamente\n"; --i; continue;
                }
                if (m <= 0) { std::cout << "Massa precisa ser > 0\n"; --i; continue; }
                masses.emplace_back(x,y,z,m);
            }
        }
        else if (mass_mode == 2) {
            // Symmetric pairs along X with given separation
            double pair_distance = prompt_double("Distância entre os pares (separação total, ex: 3.0)", 3.0);
            double mass_value = prompt_double("Valor de massa padrão para cada corpo", 5.0);
            // generate pairs until n_masses reached
            std::size_t created = 0;
            double spacing = 0.0; // offset multiplier between successive pairs
            while (created < n_masses) {
                double half = pair_distance / 2.0 + spacing;
                if (created + 1 < n_masses) {
                    masses.emplace_back(-half, 0.0, 0.0, mass_value);
                    masses.emplace_back(half, 0.0, 0.0, mass_value);
                    created += 2;
                } else {
                    // odd leftover: put central mass
                    masses.emplace_back(0.0, 0.0, 0.0, mass_value);
                    created += 1;
                }
                spacing += 0.5; // increase separation for next pair if needed
            }
        }
        else if (mass_mode == 3) {
            // Random placement with min distance
            double min_d = prompt_double("Distância mínima entre massas", 0.5);
            double mass_value = prompt_double("Valor de massa padrão", 5.0);
            std::mt19937_64 rng(std::random_device{}());
            std::uniform_real_distribution<double> ux(grid.origin()[0], grid.origin()[0] + (nx-1)*dx);
            std::uniform_real_distribution<double> uy(grid.origin()[1], grid.origin()[1] + (ny-1)*dy);
            std::uniform_real_distribution<double> uz(grid.origin()[2], grid.origin()[2] + (nz-1)*dz);

            const std::size_t max_attempts_per_point = 5000;
            for (std::size_t i = 0; i < n_masses; ++i) {
                bool placed = false;
                for (std::size_t attempt = 0; attempt < max_attempts_per_point; ++attempt) {
                    double x = ux(rng), y = uy(rng), z = uz(rng);
                    bool ok = true;
                    for (auto &t : masses) {
                        double ox2,oy2,oz2,m2;
                        std::tie(ox2,oy2,oz2,m2) = t;
                        if (dist3(x,y,z,ox2,oy2,oz2) < min_d) { ok = false; break; }
                    }
                    if (ok) { masses.emplace_back(x,y,z,mass_value); placed = true; break; }
                }
                if (!placed) {
                    std::cerr << "Falha ao posicionar massa " << i+1 << " com min_distance=" << min_d << "\n";
                    return 1;
                }
            }
        }
        else { // file
            std::string path;
            std::cout << "Caminho do arquivo (ex: masses.txt): ";
            std::getline(std::cin, path);
            if (path.empty()) { std::cerr << "Arquivo não informado\n"; return 1; }
            std::ifstream ifs(path);
            if (!ifs) { std::cerr << "Não foi possível abrir o arquivo\n"; return 1; }
            std::string line;
            while (std::getline(ifs, line) && masses.size() < n_masses) {
                if (line.empty() || line[0] == '#') continue;
                std::istringstream iss(line);
                double x,y,z,m; if (iss >> x >> y >> z >> m) {
                    if (m <= 0) continue;
                    masses.emplace_back(x,y,z,m);
                }
            }
            if (masses.empty()) { std::cerr << "Nenhuma massa lida do arquivo\n"; return 1; }
        }

        // Add masses to gravity object and report
        for (auto &t : masses) {
            double x,y,z,m; std::tie(x,y,z,m) = t;
            gravity.add_mass(x,y,z,m);
        }

        std::cout << "\nMassas adicionadas (" << gravity.mass_count() << "):\n";
        for (std::size_t i = 0; i < gravity.mass_count(); ++i) {
            const auto &ma = gravity.mass(i);
            std::cout << "  [" << i << "] (" << ma.x << ", " << ma.y << ", " << ma.z << ") m=" << ma.value << "\n";
        }
        std::cout << "  Massa total do sistema: " << gravity.total_mass() << "\n";
        if (gravity.has_overlapping_masses()) {
            std::cout << "   ⚠️  Aviso: Massas muito próximas (verifique softening)\n";
        }
        std::cout << "\n";

        // ---------------- compute ----------------
        std::cout << "4. Calculando campo gravitacional...\n";
        std::cout << "   • Iniciando cálculo sobre " << grid.total_points() << " pontos...\n";
        std::cout << "   • Processando (isso pode levar alguns segundos)...\n";
        gravity.compute(field);
        std::cout << "   • Cálculo concluído!\n";

        // stats
        const double max_mag = field.max_magnitude();
        const double min_mag = field.min_positive_magnitude();
        const double norm_l2 = field.norm_l2();

        std::cout << "   • Estatísticas do campo:\n";
        std::cout << "     Magnitude máxima: " << max_mag << "\n";
        std::cout << "     Magnitude mínima (positiva): " << min_mag << "\n";
        std::cout << "     Norma L2 do campo: " << norm_l2 << "\n\n";

        // simple invalid check
        bool has_invalid = false;
        for (std::size_t i = 0; i < 10 && i < field.size(); ++i) {
            const auto vec = field.get_vector(i);
            if (!std::isfinite(vec[0]) || !std::isfinite(vec[1]) || !std::isfinite(vec[2])) {
                has_invalid = true; break;
            }
        }
        if (has_invalid) std::cout << "   ⚠️  Campo contém NaN/inf\n";

        // ---------------- export VTI ----------------
        std::cout << "5. Exportando dados para visualização...\n";
        const std::string filename = "output.vti";
        VTIWriter::set_precision(6);
        VTIWriter::enable_magnitude(true);
        VTIWriter::enable_vector_field(true);
        VTIWriter::write(filename, grid, field);
        std::cout << "   ✓ Arquivo gerado: " << filename << "\n\n";

        // ---------------- quick validation prints ----------------
        std::cout << "6. Validação rápida (exemplos):\n";
        for (std::size_t i = 0; i < 3 && i < field.size(); ++i) {
            const auto coords = grid.linear_to_physical(i);
            const auto vec = field.get_vector(i);
            const double mag = field.magnitude(i);
            std::cout << " Ponto " << i << " @ (" << coords[0] << ", " << coords[1] << ", " << coords[2] << "):\n";
            std::cout << "   Vetor: (" << vec[0] << ", " << vec[1] << ", " << vec[2] << ")  Magnitude: " << mag << "\n";
        }
        std::cout << "\n========================================\n";
        std::cout << "   Simulação concluída com sucesso!\n";
        std::cout << "   Arquivo de saída: " << filename << "\n";
        std::cout << "========================================\n";

        return 0;
    }
    catch (const std::exception &e) {
        std::cerr << "❌ ERRO: " << e.what() << "\n";
        return 1;
    }
    catch (...) {
        std::cerr << "❌ ERRO DESCONHECIDO\n";
        return 2;
    }
}
