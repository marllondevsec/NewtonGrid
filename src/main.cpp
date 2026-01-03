#include <iostream>
#include <exception>
#include <string>

// Inclui os cabeçalhos do NewtonGrid
#include "Grid.hpp"
#include "FieldData.hpp"
#include "GravitationalField.hpp"
#include "VTIWriter.hpp"

using namespace NewtonGrid;

int main() {
    std::cout << "========================================\n";
    std::cout << "   NewtonGrid: Simulação Gravitacional 3D\n";
    std::cout << "========================================\n\n";
    
    try {
        // =====================================================================
        // 1. CONFIGURAÇÃO DO DOMÍNIO ESPACIAL
        // =====================================================================
        std::cout << "1. Criando grid cartesiano 3D...\n";
        
        // Dimensões do grid (50x50x50 pontos)
        const std::size_t nx = 50;
        const std::size_t ny = 50;
        const std::size_t nz = 50;
        
        // Espaçamento físico (unidades adimensionais)
        const double dx = 0.1;  // Espaçamento em x
        const double dy = 0.1;  // Espaçamento em y
        const double dz = 0.1;  // Espaçamento em z
        
        // Origem do sistema de coordenadas (centro do domínio)
        const double ox = -(nx * dx) / 2.0;  // -2.5
        const double oy = -(ny * dy) / 2.0;  // -2.5
        const double oz = -(nz * dz) / 2.0;  // -2.5
        
        // Cria o grid 3D
        Grid grid(nx, ny, dx, dy, ox, oy, dz, oz);
        
        std::cout << "   • Dimensões: " << nx << " × " << ny << " × " << nz << " pontos\n";
        std::cout << "   • Espaçamento: " << dx << " × " << dy << " × " << dz << "\n";
        std::cout << "   • Origem: (" << ox << ", " << oy << ", " << oz << ")\n";
        std::cout << "   • Pontos totais: " << grid.total_points() << "\n";
        std::cout << "   • Centro físico: (" << grid.center()[0] << ", "
                  << grid.center()[1] << ", " << grid.center()[2] << ")\n";
        std::cout << "   ✓ Grid criado com sucesso!\n\n";
        
        // =====================================================================
        // 2. ALOCAÇÃO DO CAMPO VETORIAL
        // =====================================================================
        std::cout << "2. Alocando campo vetorial...\n";
        
        // Cria o FieldData associado ao grid (inicializa com zeros)
        FieldData field(grid, true);
        
        std::cout << "   • FieldData alocado com " << field.size() << " pontos\n";
        std::cout << "   • Layout de memória: SoA (3 × " << field.size() << " doubles)\n";
        std::cout << "   ✓ Campo vetorial pronto para cálculo!\n\n";
        
        // =====================================================================
        // 3. CONFIGURAÇÃO DO SISTEMA FÍSICO
        // =====================================================================
        std::cout << "3. Configurando sistema gravitacional...\n";
        
        // Cria o calculador de campo gravitacional
        // G = 1.0 (unidades adimensionais)
        // Softening = 1e-3 (para estabilidade numérica)
        GravitationalField gravity(1.0, 1e-3);
        
        std::cout << "   • Constante gravitacional: " << gravity.gravitational_constant() << "\n";
        std::cout << "   • Softening: " << gravity.softening() << "\n";
        
        // Adiciona massas pontuais ao sistema
        // Massa 1: No lado negativo do eixo X
        const double mass1_value = 5.0;
        const double mass1_x = -1.5;
        const double mass1_y = 0.0;
        const double mass1_z = 0.0;
        
        // Massa 2: No lado positivo do eixo X
        const double mass2_value = 5.0;
        const double mass2_x = 1.5;
        const double mass2_y = 0.0;
        const double mass2_z = 0.0;
        
        gravity.add_mass(mass1_x, mass1_y, mass1_z, mass1_value);
        gravity.add_mass(mass2_x, mass2_y, mass2_z, mass2_value);
        
        std::cout << "   • Massas adicionadas:\n";
        std::cout << "     Massa 1: (" << mass1_x << ", " << mass1_y << ", " << mass1_z 
                  << ") com valor " << mass1_value << "\n";
        std::cout << "     Massa 2: (" << mass2_x << ", " << mass2_y << ", " << mass2_z 
                  << ") com valor " << mass2_value << "\n";
        std::cout << "   • Massa total do sistema: " << gravity.total_mass() << "\n";
        
        if (gravity.has_overlapping_masses()) {
            std::cout << "   ⚠️  Aviso: Massas estão muito próximas (possível singularidade)\n";
        }
        
        std::cout << "   ✓ Sistema gravitacional configurado!\n\n";
        
        // =====================================================================
        // 4. CÁLCULO DO CAMPO GRAVITACIONAL
        // =====================================================================
        std::cout << "4. Calculando campo gravitacional...\n";
        std::cout << "   • Iniciando cálculo sobre " << grid.total_points() << " pontos...\n";
        
        // Marca o tempo de início (conceitual - poderia usar <chrono>)
        std::cout << "   • Processando (isso pode levar alguns segundos)...\n";
        
        // Calcula o campo gravitacional
        gravity.compute(field);
        
        std::cout << "   • Cálculo concluído!\n";
        
        // Calcula algumas estatísticas do campo
        const double max_mag = field.max_magnitude();
        const double min_mag = field.min_positive_magnitude();
        const double norm_l2 = field.norm_l2();
        
        std::cout << "   • Estatísticas do campo:\n";
        std::cout << "     Magnitude máxima: " << max_mag << "\n";
        std::cout << "     Magnitude mínima (positiva): " << min_mag << "\n";
        std::cout << "     Norma L2 do campo: " << norm_l2 << "\n";
        
        // Verifica valores inválidos (NaN ou infinito)
        bool has_invalid = false;
        for (std::size_t i = 0; i < 10 && i < field.size(); ++i) {
            const auto vec = field.get_vector(i);
            if (!std::isfinite(vec[0]) || !std::isfinite(vec[1]) || !std::isfinite(vec[2])) {
                has_invalid = true;
                break;
            }
        }
        
        if (has_invalid) {
            std::cout << "   ⚠️  Aviso: Campo contém valores não-finitos (NaN/inf)\n";
        }
        
        std::cout << "   ✓ Campo gravitacional calculado com sucesso!\n\n";
        
        // =====================================================================
        // 5. EXPORTAÇÃO PARA VISUALIZAÇÃO
        // =====================================================================
        std::cout << "5. Exportando dados para visualização...\n";
        
        const std::string filename = "output.vti";
        
        // Configura o VTIWriter
        VTIWriter::set_precision(6);  // 6 dígitos são suficientes para visualização
        VTIWriter::enable_magnitude(true);
        VTIWriter::enable_vector_field(true);
        
        std::cout << "   • Criando arquivo: " << filename << "\n";
        std::cout << "   • Formato: VTK ImageData (.vti)\n";
        std::cout << "   • Compatível com: ParaView, VisIt, VTK\n";
        
        // Exporta para VTI
        VTIWriter::write(filename, grid, field);
        
        std::cout << "   ✓ Arquivo gerado com sucesso!\n\n";
        
        // =====================================================================
        // 6. INFORMAÇÕES PARA VISUALIZAÇÃO
        // =====================================================================
        std::cout << "6. Instruções para visualização:\n";
        std::cout << "   • Abra o arquivo '" << filename << "' no ParaView\n";
        std::cout << "   • Aplique o filtro 'Glyph' para visualizar vetores\n";
        std::cout << "   • Ajuste a escala do glyph para melhor visualização\n";
        std::cout << "   • Use 'Warp By Vector' para deformação por magnitude\n";
        std::cout << "   • Aplique 'Contour' para visualizar superfícies equipotenciais\n";
        std::cout << "   • O campo vetorial 'GravitationalField' contém (gx, gy, gz)\n";
        std::cout << "   • O campo escalar 'Magnitude' contém a norma do campo\n\n";
        
        // =====================================================================
        // 7. VALIDAÇÃO RÁPIDA
        // =====================================================================
        std::cout << "7. Validação rápida dos resultados:\n";
        
        // Mostra alguns valores de exemplo
        std::cout << "   • Exemplos de valores do campo (primeiros 3 pontos):\n";
        for (std::size_t i = 0; i < 3 && i < field.size(); ++i) {
            const auto coords = grid.linear_to_physical(i);
            const auto vec = field.get_vector(i);
            const double mag = field.magnitude(i);
            
            std::cout << "     Ponto " << i << " @ (" 
                      << coords[0] << ", " << coords[1] << ", " << coords[2] << "):\n";
            std::cout << "       Vetor: (" << vec[0] << ", " << vec[1] << ", " << vec[2] << ")\n";
            std::cout << "       Magnitude: " << mag << "\n";
        }
        
        // Verifica propriedades de simetria
        std::cout << "   • Verificando simetrias (pontos simétricos em relação à origem):\n";
        
        // Ponto no centro
        const auto center_indices = grid.physical_to_discrete(0.0, 0.0, 0.0);
        const std::size_t center_idx = grid.linear_index(
            center_indices[0], center_indices[1], center_indices[2]);
        
        // Ponto simétrico no eixo X
        const auto sym_indices = grid.physical_to_discrete(1.0, 0.0, 0.0);
        const std::size_t sym_idx = grid.linear_index(
            sym_indices[0], sym_indices[1], sym_indices[2]);
        
        if (center_idx < field.size() && sym_idx < field.size()) {
            const auto center_vec = field.get_vector(center_idx);
            const auto sym_vec = field.get_vector(sym_idx);
            
            std::cout << "     Centro (0,0,0): campo ≈ (" 
                      << center_vec[0] << ", " << center_vec[1] << ", " << center_vec[2] << ")\n";
            std::cout << "     Ponto (1,0,0): campo ≈ (" 
                      << sym_vec[0] << ", " << sym_vec[1] << ", " << sym_vec[2] << ")\n";
            
            // Em um sistema simétrico, o campo no centro deve ser próximo de zero
            const double center_mag = field.magnitude(center_idx);
            if (center_mag < 1e-3) {
                std::cout << "     ✓ Campo no centro é pequeno (simetria preservada)\n";
            }
        }
        
        std::cout << "\n";
        std::cout << "========================================\n";
        std::cout << "   Simulação concluída com sucesso!\n";
        std::cout << "   Arquivo de saída: " << filename << "\n";
        std::cout << "========================================\n";
        
        return 0;
        
    } catch (const std::exception& e) {
        // Tratamento de exceções
        std::cerr << "\n❌ ERRO: " << e.what() << "\n\n";
        std::cerr << "Stack trace conceitual:\n";
        std::cerr << "1. Verifique se todos os cabeçalhos estão no include path\n";
        std::cerr << "2. Confira parâmetros do grid (dimensões positivas, espaçamento > 0)\n";
        std::cerr << "3. Verifique se há espaço em disco para o arquivo de saída\n";
        std::cerr << "4. Certifique-se que FieldData foi alocado antes do cálculo\n";
        std::cerr << "5. Massas devem ter valores positivos\n";
        
        return 1;
    } catch (...) {
        std::cerr << "\n❌ ERRO DESCONHECIDO: Exceção não tratada\n";
        return 2;
    }
}
