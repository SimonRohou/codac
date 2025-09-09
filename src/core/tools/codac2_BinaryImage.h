/** 
 *  \file codac2_BinaryImage.h
 * ----------------------------------------------------------------------------
 *  \date       2025
 *  \author     Simon Rohou
 *  \copyright  Copyright 2025 Codac Team
 *  \license    GNU Lesser General Public License (LGPL)
 */

#pragma once

#include <vector>
#include <fstream>
#include "codac2_serialization.h"
#include "codac2_Grid.h"

namespace codac2
{
  using ImageDataType = uint32_t;
  // The items of the binary image are not only booleans, they can be long
  // integers resulting from the associated integral image. So they are coded
  // in uint32_t.

  template<size_t N>
  class BinaryImage : public Grid
  {
    public:

      explicit BinaryImage(const Eigen::Matrix<size_t,-1,1>& size, const Vector& dx, const Vector& origin)
        : BinaryImage(size, dx, origin, Eigen::Matrix<int,-1,1>::Zero(N))
      { }

      explicit BinaryImage(const Eigen::Matrix<size_t,-1,1>& size, const Vector& dx, const Vector& origin, const Eigen::Matrix<int,-1,1>& local_shift)
        : Grid(dx, origin, local_shift), _size(size)
      {
        assert_release(size.size() == N);
        assert_release(all_same_size(size, dx, origin, local_shift));

        _array_offsets[0] = 1;
        for(Index i = 1 ; i < (Index)N ; i++)
        {
          _array_offsets[i] = _array_offsets[i-1]*_size[i-1];
          if(i == N-1)
            _data.resize(_array_offsets[i]*_size[i], 0); // all elements are set to 0
        }
      }

      BinaryImage(const std::string& file, bool load_data = true)
        : Grid(), _file_name(file)
      {
        std::ifstream f;
        f.open(file, std::ios::binary);
        if(!f.is_open())
          throw std::runtime_error("Unable to load BinaryImage: " + file);

        Grid::deserialize(f);
        codac2::deserialize(f, _size);
        codac2::deserialize(f, _array_offsets);
        codac2::deserialize(f, _computed_as_integrale);

        assert(_size.size() == N);
        assert(all_same_size(_size, _dx, _origin, _local_shift));

        BinaryImage<N>::_data_start_pos = f.tellg();

        if(load_data)
          this->load_data(f);

        f.close();
      }

      void load_data()
      {
        assert(BinaryImage<N>::_data_start_pos > 0);

        std::ifstream f;
        f.open(_file_name, std::ios::binary);
        if(!f.is_open())
          throw std::runtime_error("Unable to load BinaryImage: " + _file_name);
        
        load_data(f);
      }

      size_t nb_voxels() const
      {
        size_t data_size = _size[0];
        for(Index i = 1 ; i < _size.size() ; i++) data_size *= _size[i];
        return data_size;
      }

      void load_data(std::ifstream& f)
      {
        f.seekg(BinaryImage<N>::_data_start_pos, f.beg);
        size_t data_size = nb_voxels();

        if(data_size != 0)
        {
          _data.resize(data_size);
          for(size_t i = 0 ; i < data_size ; i++)
            f.read((char*)&_data[i], sizeof(ImageDataType));
          assert(_data.size() == data_size);
        }
      }

      void clear_data()
      {
        _data.clear();
      }

      bool data_has_been_loaded() const
      {
        return !_data.empty() || nb_voxels() == 0;
        // nb_voxels() == 0 <=> _data.empty()
      }

      const Eigen::Matrix<size_t,N,1>& size() const
      {
        return _size;
      }

      void serialize() const
      {
        serialize(_file_name);
      }

      void serialize(const std::string& file) const
      {
        assert(file != "");
        std::ofstream f;
        f.open(file, std::ios::binary);
        if(!f.is_open())
          throw std::runtime_error("Unable to serialize: " + file);

        Grid::serialize(f);
        codac2::serialize(f, _size);
        codac2::serialize(f, _array_offsets);
        codac2::serialize(f, _computed_as_integrale);
        f.write(reinterpret_cast<const char*>(_data.data()), sizeof(ImageDataType)*_data.size());
        // todo: use a smaller type than ImageDataType for serialization
        
        f.close();
      }

      GridBox to_grid(const IntervalVector& x) const
      {
        auto g = Grid::to_grid(x);

        for(size_t i = 0 ; i < dim() ; i++)
          if(x[i].ub() == real_box()[i].ub()) // upper border of the image
            g[i] = Interval(g[i].lb(),local_indices_box()[i].ub());

        return g;
      }

      ImageDataType& operator()(const std::array<int,N>& x)
      {
        return const_cast<ImageDataType&>(static_cast<const BinaryImage<N>&>(*this).operator()(x));
      }

      const ImageDataType& operator()(const std::array<int,N>& x) const
      {
        if(*std::min_element(std::begin(x),std::end(x)) < 0)
          return _zero; // used for integral computation

        assert(local_indices_box().contains(to_vector<N>(x)));
        size_t array_index = _array_offsets[0]*x[0];
        for(size_t i = 1 ; i < N ; i++)
          array_index += _array_offsets[i]*x[i];

        assert(array_index < _data.size());
        return _data[array_index];
      }

      GridBox local_indices_box() const
      {
        GridBox b(_size.size());
        for(Index i = 0 ; i < _size.size() ; i++)
          if(_size[i] != 0)
            b[i] = { 0, _size[i]-1. };

        assert(b.has_integer_bounds());
        return b;
      }

      GridBox absolute_indices_box() const
      {
        GridBox i = GridBox(Grid::local_shift().template cast<double>()) + local_indices_box();
        assert(i.has_integer_bounds());
        return i;
      }

      IntervalVector real_box() const
      {
        return to_reals(local_indices_box());
      }

      void convert_as_integrale_img()
      {
        static_assert((N == 2 || N == 3) && "dim > 3 cases are not handled");
        // Note: for higher dimensions, one can refer to the publication
        // "A note on the computation of high-dimensional integral images",
        // Ernesto Tapia, Pattern Recognition Letters
        // https://doi.org/10.1016/j.patrec.2010.10.007

        if(_computed_as_integrale)
          return;
        _computed_as_integrale = true;
        
        if constexpr(N == 2)
        {
          for(int i = 0 ; i < (int)_size[0] ; i++)
            for(int j = 0 ; j < (int)_size[1] ; j++)
            {
              uint64_t s = (*this)({i-1,j}) + (*this)({i,j-1}) + (*this)({i,j}) - (*this)({i-1,j-1});
              assert(s < std::numeric_limits<ImageDataType>::max());
              (*this)({i,j}) = (ImageDataType)s;
            }
        }

        else if constexpr(N == 3)
        {
          for(int i = 0 ; i < (int)_size[0] ; i++)
            for(int j = 0 ; j < (int)_size[1] ; j++)
              for(int k = 0 ; k < (int)_size[2] ; k++)
              {
                uint64_t s = (*this)({i-1,j,k}) + (*this)({i,j-1,k}) + (*this)({i,j,k-1})
                  + (*this)({i-1,j-1,k-1}) + (*this)({i,j,k}) - (*this)({i-1,j,k-1})
                  - (*this)({i,j-1,k-1}) - (*this)({i-1,j-1,k});
                assert(s < std::numeric_limits<ImageDataType>::max()
                  && "integrale value exceeds the maximum value of ImageDataType");
                (*this)({i,j,k}) = (ImageDataType)s;
              }
        }
      }

      void contract_indices(GridBox& x) const
      {
        assert(x.has_integer_bounds());
        if(_data.size() == 0 || !x.is_subset(local_indices_box()) || x.is_empty())
          return; // no contraction

        assert(data_has_been_loaded());
      
        if(nb_pixels_in(x) == 0)
        {
          x.set_empty();
          return;
        }

        #if 1 // dichotomic algorithm

          for(size_t k = 0 ; k < N ; k++)
          {
            int a, b;
            auto x_ = x;

            a = x[k].lb(); b = x[k].ub();
            do {
              int i = std::floor((a+b)/2.);
              x_[k] = Interval(x[k].lb(),i);
              if(nb_pixels_in(x_) == 0)
                a = i+1;
              else
                b = i;
            } while(b!=a);
            x[k] = Interval(a,x[k].ub());

            a = x[k].lb(); b = x[k].ub();
            do {
              int i = std::ceil((a+b)/2.);
              x_[k] = Interval(i,x[k].ub());
              if(nb_pixels_in(x_) == 0)
                b = i-1;
              else
                a = i;
            } while(b!=a);
            x[k] = Interval(x[k].lb(),a);
          }

        #else // iterative algorithm

          for(size_t k = 0 ; k < N ; k++)
          {
            int imin = x[k].lb(), imax = x[k].ub(); 
            for(auto i = imin ; i <= imax ; i++)
            {
              auto x_ = x;
              x_[k] = Interval(imin,i);
              if(nb_pixels_in(x_) > 0)
              {
                imin = i;
                break;
              }
            }
            x[k] = Interval(imin,imax);

            for(auto i = imax ; i >= imin ; i--)
            {
              auto x_ = x;
              x_[k] = Interval(i,imax);
              if(nb_pixels_in(x_) > 0)
              {
                imax = i;
                break;
              }
            }
            x[k] = Interval(imin,imax);
          }

        #endif
      }

      unsigned int nb_pixels_in(const GridBox& x) const
      {
        assert(x.has_integer_bounds());
        assert(_computed_as_integrale && "conversion must be done beforehand");
        assert(!local_indices_box().is_empty());
        assert(!x.is_unbounded());

        if constexpr(N == 2)
        {
          int xmin = x[0].lb(), xmax = x[0].ub();
          int ymin = x[1].lb(), ymax = x[1].ub();

          unsigned int n1 = (*this)({xmax,   ymax  });
          unsigned int n2 = (*this)({xmax,   ymin-1});
          unsigned int n3 = (*this)({xmin-1, ymax  });
          unsigned int n4 = (*this)({xmin-1, ymin-1});
          
          return n1-n2-n3+n4;
        }

        else if constexpr(N == 3)
        {
          int xmin = x[0].lb(), xmax = x[0].ub();
          int ymin = x[1].lb(), ymax = x[1].ub();
          int zmin = x[2].lb(), zmax = x[2].ub();

          unsigned int n8 = (*this)({xmax,   ymax,   zmax  });
          unsigned int n5 = (*this)({xmin-1, ymin-1, zmax  });
          unsigned int n6 = (*this)({xmin-1, ymax,   zmax  });
          unsigned int n7 = (*this)({xmax,   ymin-1, zmax  });
          unsigned int n4 = (*this)({xmax,   ymax,   zmin-1});
          unsigned int n1 = (*this)({xmin-1, ymin-1, zmin-1});
          unsigned int n2 = (*this)({xmin-1, ymax,   zmin-1});
          unsigned int n3 = (*this)({xmax,   ymin-1, zmin-1});

          return (n8+n5-n6-n7)-(n4+n1-n2-n3);
        }

        else
          static_assert(N >= 2 && N <= 3 && "dim > 3 cases are not handled");
      }

    protected:

      std::string _file_name;
      bool _computed_as_integrale = false;
      
      Eigen::Matrix<size_t,N,1> _size, _array_offsets;
      std::vector<ImageDataType> _data;
      ImageDataType _zero = 0;

      inline static std::streampos _data_start_pos = -1;
  };
}