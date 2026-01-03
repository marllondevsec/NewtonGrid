#ifndef NEWTONGRID_FIELDDATA_HPP
#define NEWTONGRID_FIELDDATA_HPP

#include "Grid.hpp"
#include <vector>
#include <array>
#include <cstddef>
#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <limits>

namespace NewtonGrid {

class FieldData {
public:
    using ComponentBuffer = std::vector<double>;
    using Vector3d = std::array<double, 3>;
    using MagnitudeBuffer = std::vector<double>;

    FieldData() noexcept = default;

    explicit FieldData(const Grid& grid, bool zero_init = true) 
        : grid_(&grid) {
        allocate(grid, zero_init);
    }

    FieldData(const FieldData&) = default;
    FieldData& operator=(const FieldData&) = default;
    FieldData(FieldData&&) noexcept = default;
    FieldData& operator=(FieldData&&) noexcept = default;
    ~FieldData() = default;

    void allocate(const Grid& grid, bool zero_init = true) {
        const std::size_t n = grid.total_points();
        if (n == 0) {
            throw std::invalid_argument("FieldData: grid com zero pontos");
        }
        
        grid_ = &grid;
        
        if (zero_init) {
            gx_.assign(n, 0.0);
            gy_.assign(n, 0.0);
            gz_.assign(n, 0.0);
        } else {
            gx_.resize(n);
            gy_.resize(n);
            gz_.resize(n);
        }
    }

    void clear() noexcept {
        gx_.clear();
        gy_.clear();
        gz_.clear();
        grid_ = nullptr;
    }

    bool is_allocated() const noexcept { return !gx_.empty(); }
    std::size_t size() const noexcept { return gx_.size(); }
    std::size_t total_values() const noexcept { return 3 * gx_.size(); }

    const Grid& grid() const {
        check_grid();
        return *grid_;
    }

    // =========================================================================
    // ACESSO POR ÍNDICE LINEAR
    // =========================================================================

    Vector3d get_vector(std::size_t idx) const {
        check_index(idx);
        return {gx_[idx], gy_[idx], gz_[idx]};
    }

    void set_vector(std::size_t idx, const Vector3d& vec) {
        check_index(idx);
        gx_[idx] = vec[0];
        gy_[idx] = vec[1];
        gz_[idx] = vec[2];
    }

    double& gx(std::size_t idx) {
        check_index(idx);
        return gx_[idx];
    }

    double gx(std::size_t idx) const {
        check_index(idx);
        return gx_[idx];
    }

    double& gy(std::size_t idx) {
        check_index(idx);
        return gy_[idx];
    }

    double gy(std::size_t idx) const {
        check_index(idx);
        return gy_[idx];
    }

    double& gz(std::size_t idx) {
        check_index(idx);
        return gz_[idx];
    }

    double gz(std::size_t idx) const {
        check_index(idx);
        return gz_[idx];
    }

    // =========================================================================
    // ACESSO POR COORDENADAS (i, j, k)
    // =========================================================================

    Vector3d get_vector(std::size_t i, std::size_t j, std::size_t k) const {
        const std::size_t idx = grid_->linear_index(i, j, k);
        return get_vector(idx);
    }

    void set_vector(std::size_t i, std::size_t j, std::size_t k, const Vector3d& vec) {
        const std::size_t idx = grid_->linear_index(i, j, k);
        set_vector(idx, vec);
    }

    double& gx(std::size_t i, std::size_t j, std::size_t k) {
        const std::size_t idx = grid_->linear_index(i, j, k);
        return gx_[idx];
    }

    double gx(std::size_t i, std::size_t j, std::size_t k) const {
        const std::size_t idx = grid_->linear_index(i, j, k);
        return gx_[idx];
    }

    double& gy(std::size_t i, std::size_t j, std::size_t k) {
        const std::size_t idx = grid_->linear_index(i, j, k);
        return gy_[idx];
    }

    double gy(std::size_t i, std::size_t j, std::size_t k) const {
        const std::size_t idx = grid_->linear_index(i, j, k);
        return gy_[idx];
    }

    double& gz(std::size_t i, std::size_t j, std::size_t k) {
        const std::size_t idx = grid_->linear_index(i, j, k);
        return gz_[idx];
    }

    double gz(std::size_t i, std::size_t j, std::size_t k) const {
        const std::size_t idx = grid_->linear_index(i, j, k);
        return gz_[idx];
    }

    // =========================================================================
    // OPERAÇÕES EM MASSA
    // =========================================================================

    void zero() noexcept {
        std::fill(gx_.begin(), gx_.end(), 0.0);
        std::fill(gy_.begin(), gy_.end(), 0.0);
        std::fill(gz_.begin(), gz_.end(), 0.0);
    }

    void fill(double value) noexcept {
        std::fill(gx_.begin(), gx_.end(), value);
        std::fill(gy_.begin(), gy_.end(), value);
        std::fill(gz_.begin(), gz_.end(), value);
    }

    void fill(const Vector3d& vec) noexcept {
        std::fill(gx_.begin(), gx_.end(), vec[0]);
        std::fill(gy_.begin(), gy_.end(), vec[1]);
        std::fill(gz_.begin(), gz_.end(), vec[2]);
    }

    template<typename Func>
    void apply(Func&& func) {
        check_grid();
        check_allocated();
        
        const std::size_t n = size();
        for (std::size_t idx = 0; idx < n; ++idx) {
            const auto [i, j, k] = grid_->discrete_indices(idx);
            const auto [x, y, z] = grid_->linear_to_physical(idx);
            const Vector3d vec = func(i, j, k, idx, x, y, z);
            gx_[idx] = vec[0];
            gy_[idx] = vec[1];
            gz_[idx] = vec[2];
        }
    }

    void copy_from(const FieldData& other) {
        if (grid_ != other.grid_) {
            throw std::invalid_argument("FieldData::copy_from: grids diferentes");
        }
        if (size() != other.size()) {
            throw std::invalid_argument("FieldData::copy_from: tamanhos diferentes");
        }
        
        gx_ = other.gx_;
        gy_ = other.gy_;
        gz_ = other.gz_;
    }

    // =========================================================================
    // CÁLCULOS
    // =========================================================================

    double magnitude(std::size_t idx) const {
        check_index(idx);
        const double x = gx_[idx];
        const double y = gy_[idx];
        const double z = gz_[idx];
        return std::sqrt(x*x + y*y + z*z);
    }

    double magnitude(std::size_t i, std::size_t j, std::size_t k) const {
        const std::size_t idx = grid_->linear_index(i, j, k);
        return magnitude(idx);
    }

    MagnitudeBuffer compute_magnitudes() const {
        check_allocated();
        
        const std::size_t n = size();
        MagnitudeBuffer mag(n);
        
        for (std::size_t i = 0; i < n; ++i) {
            const double x = gx_[i];
            const double y = gy_[i];
            const double z = gz_[i];
            mag[i] = std::sqrt(x*x + y*y + z*z);
        }
        
        return mag;
    }

    double max_magnitude() const {
        check_allocated();
        
        double max_sq = 0.0;
        const std::size_t n = size();
        
        for (std::size_t i = 0; i < n; ++i) {
            const double x = gx_[i];
            const double y = gy_[i];
            const double z = gz_[i];
            const double mag_sq = x*x + y*y + z*z;
            if (mag_sq > max_sq) {
                max_sq = mag_sq;
            }
        }
        
        return std::sqrt(max_sq);
    }

    double min_positive_magnitude() const {
        check_allocated();
        
        double min_val = std::numeric_limits<double>::max();
        const std::size_t n = size();
        bool found = false;
        
        for (std::size_t i = 0; i < n; ++i) {
            const double x = gx_[i];
            const double y = gy_[i];
            const double z = gz_[i];
            const double mag_sq = x*x + y*y + z*z;
            
            if (mag_sq > 0.0) {
                const double mag = std::sqrt(mag_sq);
                if (mag < min_val) {
                    min_val = mag;
                    found = true;
                }
            }
        }
        
        return found ? min_val : 0.0;
    }

    double norm_l2() const {
        check_allocated();
        
        double sum = 0.0;
        const std::size_t n = size();
        
        for (std::size_t i = 0; i < n; ++i) {
            const double x = gx_[i];
            const double y = gy_[i];
            const double z = gz_[i];
            sum += x*x + y*y + z*z;
        }
        
        return std::sqrt(sum);
    }

    // =========================================================================
    // ACESSO AOS BUFFERS (para visualização)
    // =========================================================================

    ComponentBuffer& buffer_x() noexcept { return gx_; }
    const ComponentBuffer& buffer_x() const noexcept { return gx_; }
    
    ComponentBuffer& buffer_y() noexcept { return gy_; }
    const ComponentBuffer& buffer_y() const noexcept { return gy_; }
    
    ComponentBuffer& buffer_z() noexcept { return gz_; }
    const ComponentBuffer& buffer_z() const noexcept { return gz_; }
    
    double* data_x() noexcept { return gx_.data(); }
    const double* data_x() const noexcept { return gx_.data(); }
    
    double* data_y() noexcept { return gy_.data(); }
    const double* data_y() const noexcept { return gy_.data(); }
    
    double* data_z() noexcept { return gz_.data(); }
    const double* data_z() const noexcept { return gz_.data(); }
    
    std::array<double*, 3> data_all() noexcept {
        return {gx_.data(), gy_.data(), gz_.data()};
    }
    
    std::array<const double*, 3> data_all() const noexcept {
        return {gx_.data(), gy_.data(), gz_.data()};
    }

    // =========================================================================
    // CONVERSÕES PARA VISUALIZAÇÃO
    // =========================================================================

    std::vector<double> to_aos() const {
        check_allocated();
        
        const std::size_t n = size();
        std::vector<double> aos(3 * n);
        
        for (std::size_t i = 0; i < n; ++i) {
            aos[3*i]     = gx_[i];
            aos[3*i + 1] = gy_[i];
            aos[3*i + 2] = gz_[i];
        }
        
        return aos;
    }

private:
    ComponentBuffer gx_;
    ComponentBuffer gy_;
    ComponentBuffer gz_;
    const Grid* grid_ = nullptr;

    void check_index(std::size_t idx) const {
        if (idx >= gx_.size()) {
            throw std::out_of_range("FieldData: índice fora dos limites");
        }
    }

    void check_allocated() const {
        if (!is_allocated()) {
            throw std::logic_error("FieldData: não alocado");
        }
    }

    void check_grid() const {
        if (grid_ == nullptr) {
            throw std::logic_error("FieldData: sem grid associado");
        }
    }
};

} // namespace NewtonGrid

#endif // NEWTONGRID_FIELDDATA_HPP
