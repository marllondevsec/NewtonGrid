#ifndef NEWTONGRID_VTIWRITER_HPP
#define NEWTONGRID_VTIWRITER_HPP


#include "../data/Grid.hpp"       
#include "../data/FieldData.hpp"     
#include <fstream>
#include <string>
#include <sstream>
#include <iomanip>
#include <vector>
#include <array>
#include <stdexcept>
#include <cstddef>
#include <cmath>
#include <limits>

namespace NewtonGrid {

/**
 * @brief Exporta dados de simulação para o formato VTK ImageData (.vti).
 * 
 * Esta classe escreve arquivos VTI (VTK Image Data) compatíveis com ParaView,
 * representando campos vetoriais definidos sobre grades cartesianas regulares.
 * 
 * Formato gerado:
 * - VTK ImageData XML (não-binário, ASCII para simplicidade)
 * - Campo vetorial 3D (gx, gy, gz) como DataArray
 * - Magnitude do campo como DataArray escalar
 * - Compatível com grades 2D e 3D
 * - Precisão Float64 (double)
 * 
 * @note Não usa biblioteca VTK, apenas gera XML manualmente.
 * @warning Para grades muito grandes, arquivos ASCII podem ser grandes.
 *          Considere usar formatação binária para produção.
 */
class VTIWriter {
public:
    /**
     * @brief Escreve um campo vetorial em formato VTI.
     * 
     * @param filename Caminho do arquivo de saída (deve terminar em .vti)
     * @param grid Grade cartesiana que define o domínio
     * @param field Dados do campo vetorial a serem exportados
     * 
     * @throws std::invalid_argument se filename estiver vazio ou não terminar em .vti
     * @throws std::runtime_error se não for possível criar/abrir o arquivo
     * @throws std::logic_error se grid e field forem incompatíveis
     */
    static void write(const std::string& filename, 
                      const Grid& grid, 
                      const FieldData& field) {
        validate_input(filename, grid, field);
        
        std::ofstream file(filename);
        if (!file.is_open()) {
            throw std::runtime_error("Não foi possível abrir arquivo: " + filename);
        }
        
        try {
            write_header(file);
            write_image_data(file, grid);
            write_piece(file, grid, field);
            file << "</VTKFile>\n";
            
            // Verificação final
            if (!file.good()) {
                throw std::runtime_error("Erro ao escrever no arquivo: " + filename);
            }
        } catch (...) {
            file.close();
            throw;
        }
    }
    
    /**
     * @brief Configura a precisão decimal para números de ponto flutuante.
     * 
     * @param precision Número de dígitos decimais (padrão: 15)
     */
    static void set_precision(int precision) {
        if (precision < 1 || precision > 30) {
            throw std::invalid_argument("Precisão deve estar entre 1 e 30");
        }
        precision_ = precision;
    }
    
    /**
     * @brief Habilita/desabilita a exportação da magnitude.
     * 
     * @param enabled true para exportar magnitude (padrão), false para omitir
     */
    static void enable_magnitude(bool enabled) noexcept {
        export_magnitude_ = enabled;
    }
    
    /**
     * @brief Habilita/desabilita a exportação do campo vetorial.
     * 
     * @param enabled true para exportar campo vetorial (padrão), false para omitir
     */
    static void enable_vector_field(bool enabled) noexcept {
        export_vector_field_ = enabled;
    }

private:
    static inline int precision_ = 15;               ///< Precisão decimal
    static inline bool export_magnitude_ = true;     ///< Exportar magnitude
    static inline bool export_vector_field_ = true;  ///< Exportar campo vetorial
    
    /**
     * @brief Valida os parâmetros de entrada.
     */
    static void validate_input(const std::string& filename, 
                               const Grid& grid, 
                               const FieldData& field) {
        // Validação do nome do arquivo
        if (filename.empty()) {
            throw std::invalid_argument("Nome do arquivo não pode ser vazio");
        }
        
        if (filename.length() < 5 || 
            filename.substr(filename.length() - 4) != ".vti") {
            throw std::invalid_argument("Nome do arquivo deve terminar com .vti");
        }
        
        // Validação do FieldData
        if (!field.is_allocated()) {
            throw std::logic_error("FieldData não está alocado");
        }
        
        // Validação de compatibilidade entre Grid e FieldData
        if (field.size() != grid.total_points()) {
            throw std::logic_error("Grid e FieldData têm tamanhos diferentes");
        }
        
        // Validação do Grid
        const auto dims = grid.dimensions();
        if (dims[0] == 0 || dims[1] == 0 || dims[2] == 0) {
            throw std::logic_error("Grid tem dimensões inválidas (zero)");
        }
    }
    
    /**
     * @brief Escreve o cabeçalho XML do arquivo VTI.
     */
    static void write_header(std::ofstream& file) {
        file << "<?xml version=\"1.0\"?>\n";
        file << "<VTKFile type=\"ImageData\" version=\"1.0\" "
             << "byte_order=\"LittleEndian\" "
             << "header_type=\"UInt64\">\n";
    }
    
    /**
     * @brief Escreve a seção ImageData com metadados do grid.
     */
    static void write_image_data(std::ofstream& file, const Grid& grid) {
        const auto dims = grid.dimensions();
        const auto origin = grid.origin();
        const auto spacing = grid.spacing();
        
        file << "  <ImageData "
             << "WholeExtent=\"0 " << dims[0] - 1 << " "
             << "0 " << dims[1] - 1 << " "
             << "0 " << dims[2] - 1 << "\" "
             << "Origin=\""
             << std::setprecision(precision_) << origin[0] << " "
             << origin[1] << " "
             << origin[2] << "\" "
             << "Spacing=\""
             << spacing[0] << " "
             << spacing[1] << " "
             << spacing[2] << "\">\n";
    }
    
    /**
     * @brief Escreve a seção Piece com os dados do campo.
     */
    static void write_piece(std::ofstream& file, 
                            const Grid& grid, 
                            const FieldData& field) {
        const auto dims = grid.dimensions();
        
        file << "    <Piece Extent=\"0 " << dims[0] - 1 << " "
             << "0 " << dims[1] - 1 << " "
             << "0 " << dims[2] - 1 << "\">\n";
        
        file << "      <PointData>\n";
        
        if (export_vector_field_) {
            write_vector_data(file, field, "GravitationalField");
        }
        
        if (export_magnitude_) {
            write_scalar_data(file, field, "Magnitude");
        }
        
        file << "      </PointData>\n";
        file << "    </Piece>\n";
        file << "  </ImageData>\n";
    }
    
    /**
     * @brief Escreve dados vetoriais como DataArray.
     */
    static void write_vector_data(std::ofstream& file, 
                                  const FieldData& field, 
                                  const std::string& name) {
        file << "        <DataArray type=\"Float64\" "
             << "Name=\"" << name << "\" "
             << "NumberOfComponents=\"3\" "
             << "format=\"ascii\">\n";
        
        // Configura precisão
        file << std::scientific << std::setprecision(precision_);
        
        // Escreve dados no formato intercalado (AoS) para VTK
        const std::size_t n = field.size();
        for (std::size_t i = 0; i < n; ++i) {
            // VTK espera formato: vx vy vz vx vy vz ...
            file << "          " 
                 << field.gx(i) << " "
                 << field.gy(i) << " "
                 << field.gz(i);
            
            // Quebra de linha a cada 3 vetores para legibilidade
            if ((i + 1) % 3 == 0 || i == n - 1) {
                file << "\n";
            } else {
                file << " ";
            }
        }
        
        file << "        </DataArray>\n";
        
        // Restaura formatação padrão
        file << std::defaultfloat;
    }
    
    /**
     * @brief Escreve dados escalares como DataArray.
     */
    static void write_scalar_data(std::ofstream& file, 
                                  const FieldData& field, 
                                  const std::string& name) {
        file << "        <DataArray type=\"Float64\" "
             << "Name=\"" << name << "\" "
             << "NumberOfComponents=\"1\" "
             << "format=\"ascii\">\n";
        
        // Configura precisão
        file << std::scientific << std::setprecision(precision_);
        
        // Calcula magnitudes
        auto magnitudes = field.compute_magnitudes();
        
        // Escreve dados escalares
        const std::size_t n = magnitudes.size();
        for (std::size_t i = 0; i < n; ++i) {
            file << "          " << magnitudes[i];
            
            // Quebra de linha a cada 6 valores para legibilidade
            if ((i + 1) % 6 == 0 || i == n - 1) {
                file << "\n";
            } else {
                file << " ";
            }
        }
        
        file << "        </DataArray>\n";
        
        // Restaura formatação padrão
        file << std::defaultfloat;
    }
    
    /**
     * @brief Verifica se um valor double é finito (não NaN ou infinito).
     */
    static bool is_finite(double value) noexcept {
        return std::isfinite(value);
    }
    
    /**
     * @brief Converte um vetor para string formatada para depuração.
     */
    static std::string vector_to_string(const std::array<double, 3>& vec) {
        std::ostringstream oss;
        oss << std::scientific << std::setprecision(precision_)
            << "(" << vec[0] << ", " << vec[1] << ", " << vec[2] << ")";
        return oss.str();
    }
};

} // namespace NewtonGrid

#endif // NEWTONGRID_VTIWRITER_HPP
