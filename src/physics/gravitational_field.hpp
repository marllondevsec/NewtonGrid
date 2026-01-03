#ifndef NEWTONGRID_GRAVITATIONALFIELD_HPP
#define NEWTONGRID_GRAVITATIONALFIELD_HPP

#include "Grid.hpp"
#include "FieldData.hpp"
#include <vector>
#include <array>
#include <cmath>
#include <stdexcept>
#include <cstddef>
#include <limits>

namespace NewtonGrid {

/**
 * @brief Representa uma massa pontual no espaço.
 * 
 * Usada como fonte para o cálculo do campo gravitacional.
 */
struct Mass {
    double x;      ///< Coordenada X da massa
    double y;      ///< Coordenada Y da massa
    double z;      ///< Coordenada Z da massa
    double value;  ///< Valor da massa (sempre positivo)
    
    /**
     * @brief Constrói uma massa pontual.
     * 
     * @param x Coordenada X
     * @param y Coordenada Y
     * @param z Coordenada Z
     * @param value Valor da massa (deve ser > 0)
     * 
     * @throws std::invalid_argument se value <= 0
     */
    constexpr Mass(double x, double y, double z, double value)
        : x(x), y(y), z(z), value(value) {
        if (value <= 0.0) {
            throw std::invalid_argument("Massa deve ser positiva");
        }
    }
};

/**
 * @brief Calcula campos gravitacionais newtonianos em grades cartesianas.
 * 
 * Esta classe implementa o cálculo do campo gravitacional gerado por
 * massas pontuais usando a lei da gravitação universal de Newton:
 * 
 *     g⃗ = -G * m * r⃗ / |r|³
 * 
 * Características:
 * - Suporta múltiplas massas pontuais
 * - Inclui softening para estabilidade numérica
 * - Cálculo determinístico e thread-safe (para leitura)
 * - Preparado para otimizações futuras (SIMD, OpenMP)
 * 
 * @note A classe não armazena o campo calculado, apenas escreve em um FieldData.
 */
class GravitationalField {
public:
    using Vector3d = std::array<double, 3>;
    
    /**
     * @brief Constrói um calculador de campo gravitacional.
     * 
     * @param G Constante gravitacional (padrão: 1.0 para unidades adimensionais)
     * @param softening Parâmetro de softening para evitar singularidades (padrão: 1e-10)
     * 
     * @throws std::invalid_argument se G <= 0 ou softening < 0
     */
    explicit GravitationalField(double G = 1.0, double softening = 1e-10)
        : G_(G), softening_squared_(softening * softening) {
        if (G <= 0.0) {
            throw std::invalid_argument("Constante gravitacional deve ser > 0");
        }
        if (softening < 0.0) {
            throw std::invalid_argument("Softening não pode ser negativo");
        }
    }
    
    // =========================================================================
    // CONFIGURAÇÃO
    // =========================================================================
    
    /**
     * @brief Define a constante gravitacional.
     * 
     * @param G Nova constante gravitacional (> 0)
     * 
     * @throws std::invalid_argument se G <= 0
     */
    void set_gravitational_constant(double G) {
        if (G <= 0.0) {
            throw std::invalid_argument("Constante gravitacional deve ser > 0");
        }
        G_ = G;
    }
    
    /**
     * @brief Retorna a constante gravitacional atual.
     */
    double gravitational_constant() const noexcept { return G_; }
    
    /**
     * @brief Define o parâmetro de softening.
     * 
     * O softening evita divisão por zero quando o ponto de cálculo está
     * muito próximo de uma massa. A fórmula modificada é:
     * 
     *     g⃗ = -G * m * r⃗ / (|r|² + ε²)^{3/2}
     * 
     * @param softening Novo valor de softening (≥ 0)
     */
    void set_softening(double softening) noexcept {
        softening_squared_ = softening * softening;
    }
    
    /**
     * @brief Retorna o parâmetro de softening atual.
     */
    double softening() const noexcept { return std::sqrt(softening_squared_); }
    
    // =========================================================================
    // GERENCIAMENTO DE MASSAS
    // =========================================================================
    
    /**
     * @brief Adiciona uma massa pontual ao sistema.
     * 
     * @param mass Objeto Mass contendo posição e valor
     */
    void add_mass(const Mass& mass) {
        masses_.push_back(mass);
    }
    
    /**
     * @brief Adiciona uma massa pontual ao sistema.
     * 
     * @param x Coordenada X
     * @param y Coordenada Y
     * @param z Coordenada Z
     * @param value Valor da massa (> 0)
     * 
     * @throws std::invalid_argument se value <= 0
     */
    void add_mass(double x, double y, double z, double value) {
        masses_.emplace_back(x, y, z, value);
    }
    
    /**
     * @brief Remove todas as massas do sistema.
     */
    void clear_masses() noexcept {
        masses_.clear();
    }
    
    /**
     * @brief Retorna o número de massas no sistema.
     */
    std::size_t mass_count() const noexcept {
        return masses_.size();
    }
    
    /**
     * @brief Retorna referência à massa na posição especificada.
     * 
     * @param index Índice da massa [0, mass_count()-1]
     * @return Referência constante à massa
     * 
     * @throws std::out_of_range se index >= mass_count()
     */
    const Mass& mass(std::size_t index) const {
        if (index >= masses_.size()) {
            throw std::out_of_range("Índice de massa fora dos limites");
        }
        return masses_[index];
    }
    
    /**
     * @brief Remove a massa na posição especificada.
     * 
     * @param index Índice da massa a ser removida
     * 
     * @throws std::out_of_range se index >= mass_count()
     */
    void remove_mass(std::size_t index) {
        if (index >= masses_.size()) {
            throw std::out_of_range("Índice de massa fora dos limites");
        }
        masses_.erase(masses_.begin() + index);
    }
    
    // =========================================================================
    // CÁLCULO DO CAMPO
    // =========================================================================
    
    /**
     * @brief Calcula o campo gravitacional em todos os pontos do grid.
     * 
     * Para cada ponto do grid, calcula a contribuição de todas as massas
     * e escreve o campo resultante no FieldData.
     * 
     * Fórmula para cada massa:
     *     g⃗ += -G * m * r⃗ / (|r|² + ε²)^{3/2}
     * 
     * @param field FieldData onde o resultado será escrito (deve estar alocado)
     * 
     * @throws std::invalid_argument se field não estiver alocado
     * @throws std::logic_error se não houver massas definidas
     */
    void compute(FieldData& field) const {
        // Validação básica
        if (!field.is_allocated()) {
            throw std::invalid_argument("FieldData não está alocado");
        }
        
        if (masses_.empty()) {
            throw std::logic_error("Nenhuma massa definida para cálculo do campo");
        }
        
        const Grid& grid = field.grid();
        const std::size_t n_points = field.size();
        
        // Zera o campo antes de calcular
        field.zero();
        
        // Para cada ponto do grid
        for (std::size_t idx = 0; idx < n_points; ++idx) {
            // Obtém coordenadas físicas do ponto
            const auto [x, y, z] = grid.linear_to_physical(idx);
            
            // Acumula contribuições de todas as massas
            Vector3d total_field = {0.0, 0.0, 0.0};
            
            for (const auto& mass : masses_) {
                // Vetor deslocamento do ponto à massa
                const double dx = x - mass.x;
                const double dy = y - mass.y;
                const double dz = z - mass.z;
                
                // Distância ao quadrado com softening
                const double r_sq = dx*dx + dy*dy + dz*dz + softening_squared_;
                
                // Evita problemas numéricos para distâncias muito grandes
                if (r_sq < std::numeric_limits<double>::min()) {
                    continue;  // Ponto coincide exatamente com a massa (após softening)
                }
                
                // Calcula fator comum: -G * m / (r_sq * sqrt(r_sq))
                const double factor = -G_ * mass.value / (r_sq * std::sqrt(r_sq));
                
                // Acumula contribuição
                total_field[0] += factor * dx;
                total_field[1] += factor * dy;
                total_field[2] += factor * dz;
            }
            
            // Armazena o campo resultante no FieldData
            field.set_vector(idx, total_field);
        }
    }
    
    /**
     * @brief Calcula o campo gravitacional usando iteração otimizada.
     * 
     * Versão alternativa que usa o método for_each_point do Grid,
     * que pode ser mais eficiente para alguns casos de uso.
     * 
     * @param field FieldData onde o resultado será escrito
     */
    void compute_optimized(FieldData& field) const {
        // Validação básica
        if (!field.is_allocated()) {
            throw std::invalid_argument("FieldData não está alocado");
        }
        
        if (masses_.empty()) {
            throw std::logic_error("Nenhuma massa definida para cálculo do campo");
        }
        
        const Grid& grid = field.grid();
        
        // Zera o campo antes de calcular
        field.zero();
        
        // Itera por todos os pontos do grid
        grid.for_each_point([&](std::size_t i, std::size_t j, std::size_t k,
                               std::size_t idx, double x, double y, double z) {
            // Acumula contribuições de todas as massas
            Vector3d total_field = {0.0, 0.0, 0.0};
            
            for (const auto& mass : masses_) {
                const double dx = x - mass.x;
                const double dy = y - mass.y;
                const double dz = z - mass.z;
                
                const double r_sq = dx*dx + dy*dy + dz*dz + softening_squared_;
                
                if (r_sq < std::numeric_limits<double>::min()) {
                    continue;
                }
                
                const double factor = -G_ * mass.value / (r_sq * std::sqrt(r_sq));
                
                total_field[0] += factor * dx;
                total_field[1] += factor * dy;
                total_field[2] += factor * dz;
            }
            
            // Armazena o campo resultante
            field.set_vector(idx, total_field);
        });
    }
    
    /**
     * @brief Calcula apenas a magnitude do campo gravitacional.
     * 
     * Versão especializada que calcula apenas a magnitude (norma) do campo
     * em cada ponto, ignorando a direção.
     * 
     * @param field FieldData onde o resultado será escrito
     * 
     * @note Esta versão é mais eficiente quando apenas a magnitude é necessária.
     */
    void compute_magnitude_only(FieldData& field) const {
        // Validação básica
        if (!field.is_allocated()) {
            throw std::invalid_argument("FieldData não está alocado");
        }
        
        if (masses_.empty()) {
            throw std::logic_error("Nenhuma massa definida para cálculo do campo");
        }
        
        const Grid& grid = field.grid();
        const std::size_t n_points = field.size();
        
        // Para cada ponto do grid
        for (std::size_t idx = 0; idx < n_points; ++idx) {
            const auto [x, y, z] = grid.linear_to_physical(idx);
            
            double total_magnitude = 0.0;
            
            // Para cada massa, calcula contribuição à magnitude
            // Nota: Esta aproximação não é exata para múltiplas massas,
            // mas é válida quando as massas estão bem separadas
            for (const auto& mass : masses_) {
                const double dx = x - mass.x;
                const double dy = y - mass.y;
                const double dz = z - mass.z;
                
                const double r_sq = dx*dx + dy*dy + dz*dz + softening_squared_;
                
                if (r_sq < std::numeric_limits<double>::min()) {
                    continue;
                }
                
                // Magnitude do campo de uma massa pontual
                total_magnitude += G_ * mass.value / r_sq;
            }
            
            // Armazena magnitude como se fosse componente X
            // (FieldData ainda é vetorial, mas usamos só X)
            field.gx(idx) = total_magnitude;
            field.gy(idx) = 0.0;
            field.gz(idx) = 0.0;
        }
    }
    
    /**
     * @brief Calcula o potencial gravitacional em vez do campo.
     * 
     * O potencial gravitacional Φ é definido como:
     *     Φ = -G * m / (|r| + ε)
     * 
     * @param field FieldData onde o resultado será escrito
     * 
     * @note O FieldData será usado para armazenar o potencial escalar
     *       na componente X, com Y e Z zeradas.
     */
    void compute_potential(FieldData& field) const {
        // Validação básica
        if (!field.is_allocated()) {
            throw std::invalid_argument("FieldData não está alocado");
        }
        
        if (masses_.empty()) {
            throw std::logic_error("Nenhuma massa definida para cálculo do potencial");
        }
        
        const Grid& grid = field.grid();
        const std::size_t n_points = field.size();
        
        // Para cada ponto do grid
        for (std::size_t idx = 0; idx < n_points; ++idx) {
            const auto [x, y, z] = grid.linear_to_physical(idx);
            
            double total_potential = 0.0;
            
            for (const auto& mass : masses_) {
                const double dx = x - mass.x;
                const double dy = y - mass.y;
                const double dz = z - mass.z;
                
                const double r = std::sqrt(dx*dx + dy*dy + dz*dz + softening_squared_);
                
                if (r < std::numeric_limits<double>::min()) {
                    continue;
                }
                
                total_potential += -G_ * mass.value / r;
            }
            
            // Armazena potencial como se fosse componente X
            field.gx(idx) = total_potential;
            field.gy(idx) = 0.0;
            field.gz(idx) = 0.0;
        }
    }
    
    // =========================================================================
    // UTILITÁRIOS
    // =========================================================================
    
    /**
     * @brief Verifica se há sobreposição significativa entre massas.
     * 
     * @param threshold Distância mínima considerada como sobreposição
     * @return true se alguma massa está mais próxima que threshold de outra
     */
    bool has_overlapping_masses(double threshold = 1e-6) const {
        const double threshold_sq = threshold * threshold;
        
        for (std::size_t i = 0; i < masses_.size(); ++i) {
            for (std::size_t j = i + 1; j < masses_.size(); ++j) {
                const auto& m1 = masses_[i];
                const auto& m2 = masses_[j];
                
                const double dx = m1.x - m2.x;
                const double dy = m1.y - m2.y;
                const double dz = m1.z - m2.z;
                
                if (dx*dx + dy*dy + dz*dz < threshold_sq) {
                    return true;
                }
            }
        }
        
        return false;
    }
    
    /**
     * @brief Calcula a massa total do sistema.
     */
    double total_mass() const noexcept {
        double total = 0.0;
        for (const auto& mass : masses_) {
            total += mass.value;
        }
        return total;
    }
    
private:
    double G_;                    ///< Constante gravitacional
    double softening_squared_;    ///< Softening ao quadrado (ε²)
    std::vector<Mass> masses_;    ///< Lista de massas pontuais
};

} // namespace NewtonGrid

#endif // NEWTONGRID_GRAVITATIONALFIELD_HPP
