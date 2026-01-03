#ifndef NEWTONGRID_GRID_HPP
#define NEWTONGRID_GRID_HPP

#include <vector>
#include <array>
#include <cstddef>
#include <stdexcept>
#include <cassert>
#include <type_traits>
#include <algorithm>
#include <cmath>

namespace NewtonGrid {

/**
 * @brief Representa um domínio cartesiano regular para simulações numéricas.
 * 
 * O Grid fornece uma estrutura de dados otimizada para cálculos científicos
 * em domínios 2D/3D, com indexação linear cache-friendly e conversões
 * eficientes entre coordenadas físicas e índices discretos.
 * 
 * Características principais:
 * - Representação interna sempre 3D (2D: nz=1, dz=0)
 * - Memória contígua e ordenamento F-contiguous (x mais rápido)
 * - Interface imutável após construção
 * - Conversões eficientes entre (i,j,k) <-> índice linear <-> (x,y,z)
 * - Compatível com VTK ImageData e formatos científicos
 */
class Grid {
public:
    // Tipos fundamentais
    using Index = std::array<std::size_t, 3>;    ///< Índices discretos (i, j, k)
    using Point = std::array<double, 3>;         ///< Coordenada física (x, y, z)
    using Spacing = std::array<double, 3>;       ///< Espaçamento físico (dx, dy, dz)
    using Extent = std::array<double, 6>;        ///< Extremos físicos [xmin, xmax, ymin, ymax, zmin, zmax]

    /**
     * @brief Constrói um grid cartesiano regular.
     * 
     * @param dimensions Número de pontos em cada direção [nx, ny, nz]
     * @param spacing Espaçamento físico entre pontos [dx, dy, dz] (deve ser > 0)
     * @param origin Coordenada física do ponto (0,0,0) [ox, oy, oz]
     * 
     * @throws std::invalid_argument Se dimensions contiver zero ou spacing contiver valor ≤ 0
     */
    Grid(const std::array<std::size_t, 3>& dimensions,
         const std::array<double, 3>& spacing,
         const std::array<double, 3>& origin)
        : dimensions_(dimensions),
          spacing_(spacing),
          origin_(origin),
          total_points_(dimensions[0] * dimensions[1] * dimensions[2]) {
        
        // Validação rigorosa de parâmetros
        validateConstruction();
        computeDerivedProperties();
    }

    /**
     * @brief Constrói um grid 2D (nz = 1, dz arbitrário mas tipicamente 1.0).
     * 
     * @param nx Número de pontos em X (≥ 1)
     * @param ny Número de pontos em Y (≥ 1)
     * @param dx Espaçamento em X (> 0)
     * @param dy Espaçamento em Y (> 0)
     * @param ox Origem X
     * @param oy Origem Y
     * @param dz Espaçamento em Z (tipicamente 1.0, > 0)
     * @param oz Origem Z (tipicamente 0.0)
     */
    Grid(std::size_t nx, std::size_t ny,
         double dx, double dy,
         double ox, double oy,
         double dz = 1.0, double oz = 0.0)
        : Grid({nx, ny, 1}, {dx, dy, dz}, {ox, oy, oz}) {}

    // Grid é imutável após construção - não permitir cópias acidentais
    Grid(const Grid&) = default;
    Grid& operator=(const Grid&) = delete;
    Grid(Grid&&) = default;
    Grid& operator=(Grid&&) = delete;

    // =========================================================================
    // GETTERS BÁSICOS (informação do domínio)
    // =========================================================================

    /** @brief Retorna as dimensões do grid [nx, ny, nz] */
    const std::array<std::size_t, 3>& dimensions() const noexcept { return dimensions_; }

    /** @brief Número de pontos na direção X */
    std::size_t nx() const noexcept { return dimensions_[0]; }
    
    /** @brief Número de pontos na direção Y */
    std::size_t ny() const noexcept { return dimensions_[1]; }
    
    /** @brief Número de pontos na direção Z */
    std::size_t nz() const noexcept { return dimensions_[2]; }

    /** @brief Espaçamento físico entre pontos [dx, dy, dz] */
    const std::array<double, 3>& spacing() const noexcept { return spacing_; }
    
    /** @brief Espaçamento em X */
    double dx() const noexcept { return spacing_[0]; }
    
    /** @brief Espaçamento em Y */
    double dy() const noexcept { return spacing_[1]; }
    
    /** @brief Espaçamento em Z */
    double dz() const noexcept { return spacing_[2]; }

    /** @brief Origem física do grid [ox, oy, oz] */
    const std::array<double, 3>& origin() const noexcept { return origin_; }
    
    /** @brief Coordenada X da origem */
    double ox() const noexcept { return origin_[0]; }
    
    /** @brief Coordenada Y da origem */
    double oy() const noexcept { return origin_[1]; }
    
    /** @brief Coordenada Z da origem */
    double oz() const noexcept { return origin_[2]; }

    /** @brief Número total de pontos no grid (nx * ny * nz) */
    std::size_t total_points() const noexcept { return total_points_; }

    /** @brief Retorna true se o grid for 2D (nz == 1) */
    bool is_2d() const noexcept { return dimensions_[2] == 1; }

    /** @brief Retorna true se o grid for 3D (nz > 1) */
    bool is_3d() const noexcept { return dimensions_[2] > 1; }

    // =========================================================================
    // GEOMETRIA DO DOMÍNIO
    // =========================================================================

    /**
     * @brief Retorna os extremos físicos do grid.
     * 
     * @return Array [xmin, xmax, ymin, ymax, zmin, zmax]
     *         onde: xmin = ox, xmax = ox + (nx-1)*dx
     */
    Extent physical_extent() const noexcept {
        return {
            origin_[0],                                    // xmin
            origin_[0] + (dimensions_[0] - 1) * spacing_[0], // xmax
            origin_[1],                                    // ymin
            origin_[1] + (dimensions_[1] - 1) * spacing_[1], // ymax
            origin_[2],                                    // zmin
            origin_[2] + (dimensions_[2] - 1) * spacing_[2]  // zmax
        };
    }

    /** @brief Centro físico do grid */
    Point center() const noexcept {
        const auto& ext = physical_extent();
        return {
            (ext[0] + ext[1]) * 0.5,
            (ext[2] + ext[3]) * 0.5,
            (ext[4] + ext[5]) * 0.5
        };
    }

    /** @brief Dimensões físicas (comprimento, largura, altura) */
    std::array<double, 3> physical_size() const noexcept {
        return {
            (dimensions_[0] - 1) * spacing_[0],
            (dimensions_[1] - 1) * spacing_[1],
            (dimensions_[2] - 1) * spacing_[2]
        };
    }

    // =========================================================================
    // CONVERSÕES DE COORDENADAS (core da classe)
    // =========================================================================

    /**
     * @brief Converte índices discretos (i,j,k) para índice linear.
     * 
     * Ordenamento F-contiguous (x mais rápido): idx = i + j*nx + k*nx*ny
     * 
     * @param i Índice em X [0, nx-1]
     * @param j Índice em Y [0, ny-1]
     * @param k Índice em Z [0, nz-1]
     * @return Índice linear no array 1D
     * 
     * @throws std::out_of_range Se algum índice estiver fora dos limites
     */
    std::size_t linear_index(std::size_t i, std::size_t j, std::size_t k) const {
        validate_indices(i, j, k);
        return i + j * dimensions_[0] + k * dimensions_[0] * dimensions_[1];
    }

    /**
     * @brief Converte índices discretos (i,j,k) para índice linear (sem validação).
     * 
     * Versão mais rápida para uso interno quando os índices já são conhecidos como válidos.
     */
    std::size_t linear_index_unchecked(std::size_t i, std::size_t j, std::size_t k) const noexcept {
        return i + j * dimensions_[0] + k * dimensions_[0] * dimensions_[1];
    }

    /**
     * @brief Converte índice linear para índices discretos (i,j,k).
     * 
     * @param idx Índice linear [0, total_points-1]
     * @return Array {i, j, k}
     * 
     * @throws std::out_of_range Se idx ≥ total_points
     */
    Index discrete_indices(std::size_t idx) const {
        if (idx >= total_points_) {
            throw std::out_of_range("Índice linear fora dos limites");
        }
        
        const std::size_t nxny = dimensions_[0] * dimensions_[1];
        const std::size_t k = idx / nxny;
        const std::size_t remainder = idx % nxny;
        const std::size_t j = remainder / dimensions_[0];
        const std::size_t i = remainder % dimensions_[0];
        
        return {i, j, k};
    }

    /**
     * @brief Converte coordenada física para índices discretos mais próximos.
     * 
     * @param x Coordenada X física
     * @param y Coordenada Y física
     * @param z Coordenada Z física
     * @return Índices discretos {i, j, k} do ponto do grid mais próximo
     */
    Index physical_to_discrete(double x, double y, double z) const noexcept {
        return {
            static_cast<std::size_t>(std::round((x - origin_[0]) / spacing_[0])),
            static_cast<std::size_t>(std::round((y - origin_[1]) / spacing_[1])),
            static_cast<std::size_t>(std::round((z - origin_[2]) / spacing_[2]))
        };
    }

    /**
     * @brief Converte índices discretos para coordenada física.
     * 
     * @param i Índice em X
     * @param j Índice em Y
     * @param k Índice em Z
     * @return Coordenada física {x, y, z} do ponto do grid
     */
    Point discrete_to_physical(std::size_t i, std::size_t j, std::size_t k) const {
        validate_indices(i, j, k);
        return {
            origin_[0] + i * spacing_[0],
            origin_[1] + j * spacing_[1],
            origin_[2] + k * spacing_[2]
        };
    }

    /**
     * @brief Converte índice linear para coordenada física.
     * 
     * @param idx Índice linear
     * @return Coordenada física {x, y, z}
     */
    Point linear_to_physical(std::size_t idx) const {
        const auto indices = discrete_indices(idx);
        return discrete_to_physical(indices[0], indices[1], indices[2]);
    }

    // =========================================================================
    // VALIDAÇÕES E VERIFICAÇÕES
    // =========================================================================

    /**
     * @brief Verifica se índices discretos estão dentro dos limites do grid.
     */
    bool contains_index(std::size_t i, std::size_t j, std::size_t k) const noexcept {
        return (i < dimensions_[0]) && (j < dimensions_[1]) && (k < dimensions_[2]);
    }

    /**
     * @brief Verifica se uma coordenada física está dentro dos limites do grid.
     * 
     * @param point Coordenada a verificar
     * @param tolerance Tolerância em unidades físicas (para pontos na borda)
     * @return true se o ponto estiver dentro do domínio
     */
    bool contains_physical_point(const Point& point, double tolerance = 1e-12) const noexcept {
        const auto [xmin, xmax, ymin, ymax, zmin, zmax] = physical_extent();
        
        return (point[0] >= xmin - tolerance) && (point[0] <= xmax + tolerance) &&
               (point[1] >= ymin - tolerance) && (point[1] <= ymax + tolerance) &&
               (point[2] >= zmin - tolerance) && (point[2] <= zmax + tolerance);
    }

    // =========================================================================
    // ITERAÇÃO E UTILITÁRIOS
    // =========================================================================

    /**
     * @brief Aplica uma função a cada ponto do grid.
     * 
     * @tparam Func Tipo da função (deve aceitar (i, j, k, idx, x, y, z))
     * @param func Função a aplicar
     */
    template<typename Func>
    void for_each_point(Func&& func) const {
        for (std::size_t k = 0; k < dimensions_[2]; ++k) {
            const double z = origin_[2] + k * spacing_[2];
            for (std::size_t j = 0; j < dimensions_[1]; ++j) {
                const double y = origin_[1] + j * spacing_[1];
                for (std::size_t i = 0; i < dimensions_[0]; ++i) {
                    const double x = origin_[0] + i * spacing_[0];
                    const std::size_t idx = linear_index_unchecked(i, j, k);
                    func(i, j, k, idx, x, y, z);
                }
            }
        }
    }

    /**
     * @brief Retorna um vetor com as coordenadas de todos os pontos do grid.
     * 
     * Útil para inicialização rápida ou exportação. Cuidado com grids grandes!
     */
    std::vector<Point> all_coordinates() const {
        std::vector<Point> coords;
        coords.reserve(total_points_);
        
        for_each_point([&coords](std::size_t, std::size_t, std::size_t,
                                 std::size_t, double x, double y, double z) {
            coords.push_back({x, y, z});
        });
        
        return coords;
    }

    /**
     * @brief Cria um grid com as mesmas dimensões mas espaçamento/origem diferente.
     * 
     * Útil para reamostragem ou mudança de sistema de coordenadas.
     */
    Grid with_spacing_and_origin(const Spacing& new_spacing, const Point& new_origin) const {
        return Grid(dimensions_, new_spacing, new_origin);
    }

private:
    // =========================================================================
    // DADOS MEMBRO (imutáveis após construção)
    // =========================================================================
    std::array<std::size_t, 3> dimensions_;  ///< [nx, ny, nz]
    std::array<double, 3> spacing_;          ///< [dx, dy, dz] (sempre > 0)
    std::array<double, 3> origin_;           ///< [ox, oy, oz]
    std::size_t total_points_;               ///< nx * ny * nz
    
    // Propriedades derivadas (cacheadas para performance)
    std::array<double, 3> inv_spacing_;      ///< [1/dx, 1/dy, 1/dz]

    // =========================================================================
    // MÉTODOS PRIVADOS
    // =========================================================================

    /**
     * @brief Valida parâmetros de construção.
     * 
     * @throws std::invalid_argument Para parâmetros inválidos
     */
    void validateConstruction() const {
        // Valida dimensões
        for (std::size_t i = 0; i < 3; ++i) {
            if (dimensions_[i] == 0) {
                throw std::invalid_argument("Dimensões do grid devem ser ≥ 1");
            }
        }
        
        // Valida espaçamento
        for (std::size_t i = 0; i < 3; ++i) {
            if (spacing_[i] <= 0.0) {
                throw std::invalid_argument("Espaçamento deve ser > 0");
            }
        }
        
        // Verificação adicional para grid 2D
        if (dimensions_[2] == 1 && std::abs(spacing_[2]) < 1e-12) {
            throw std::invalid_argument("Grid 2D deve ter dz > 0 (tipicamente 1.0)");
        }
    }

    /**
     * @brief Computa e cacheia propriedades derivadas.
     */
    void computeDerivedProperties() noexcept {
        for (std::size_t i = 0; i < 3; ++i) {
            inv_spacing_[i] = 1.0 / spacing_[i];
        }
    }

    /**
     * @brief Valida índices discretos.
     * 
     * @throws std::out_of_range Se algum índice estiver fora dos limites
     */
    void validate_indices(std::size_t i, std::size_t j, std::size_t k) const {
        if (i >= dimensions_[0] || j >= dimensions_[1] || k >= dimensions_[2]) {
            throw std::out_of_range("Índices discretos fora dos limites do grid");
        }
    }
};

// =============================================================================
// FUNÇÕES AUXILIARES (não-membro)
// =============================================================================

/**
 * @brief Compara dois grids por igualdade (dentro de tolerância).
 */
inline bool grids_equal(const Grid& a, const Grid& b, double tolerance = 1e-12) noexcept {
    const auto& dim_a = a.dimensions();
    const auto& dim_b = b.dimensions();
    
    // Compara dimensões
    for (int i = 0; i < 3; ++i) {
        if (dim_a[i] != dim_b[i]) return false;
    }
    
    // Compara espaçamento e origem
    const auto& sp_a = a.spacing();
    const auto& sp_b = b.spacing();
    const auto& org_a = a.origin();
    const auto& org_b = b.origin();
    
    for (int i = 0; i < 3; ++i) {
        if (std::abs(sp_a[i] - sp_b[i]) > tolerance) return false;
        if (std::abs(org_a[i] - org_b[i]) > tolerance) return false;
    }
    
    return true;
}

/**
 * @brief Cria um grid centrado na origem.
 * 
 * @param dimensions Dimensões [nx, ny, nz]
 * @param spacing Espaçamento [dx, dy, dz]
 * @return Grid centrado em (0,0,0)
 */
inline Grid centered_grid(const std::array<std::size_t, 3>& dimensions,
                         const std::array<double, 3>& spacing) {
    const double ox = -0.5 * (dimensions[0] - 1) * spacing[0];
    const double oy = -0.5 * (dimensions[1] - 1) * spacing[1];
    const double oz = -0.5 * (dimensions[2] - 1) * spacing[2];
    
    return Grid(dimensions, spacing, {ox, oy, oz});
}

} // namespace NewtonGrid

#endif // NEWTONGRID_GRID_HPP
