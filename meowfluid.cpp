#include <cstdio>
#include <assert.h>
#include <string>
#include <iostream>
#include <array>
#include <algorithm>
#include <array>
#include <random>
#include <variant>
#include <vector>
#include <unordered_map>
#include <memory>
#include <fstream>
#include <type_traits>
#include <cstdint>
#include <signal.h>
#include <stdio.h>
#include <unistd.h>


using namespace std;

constexpr std::array<pair<int, int>, 4> deltas{{{-1, 0}, {1, 0}, {0, -1}, {0, 1}}};


struct InitArgs {
  double G, RHO_F, RHO_A;
  std::vector<std::vector<char>> tmp_field;
} initArgs;

template<size_t N, bool fast>
struct FixedDataType;

template<size_t N>
struct FixedDataType<N, true> {
    using type = std::conditional_t<
        (N <= 8),
        int_fast8_t,
        std::conditional_t<
            (N <= 16),
            int_fast16_t,

#           ifdef INT32_MAX
            std::conditional_t<
                (N <= 32),
                int_fast32_t,
#               ifdef INT64_MAX
                std::conditional_t<
                    (N <= 64),
                    int_fast64_t,
#                   ifdef INT128_MAX
                    int_fast128_t
#                   else
                    double
#                   endif
                >
#               else
                double
#               endif
            >
#        else
           double
#        endif
        >
    >;
};
#ifdef INT8_MAX
template<>
struct FixedDataType<8, false> {
  using type = int8_t;
};
#endif
#ifdef INT16_MAX
template<>
struct FixedDataType<16, false> {
  using type = int16_t;
};
#endif
#ifdef INT32_MAX
template<>
struct FixedDataType<32, false> {
  using type = int32_t;
};
#endif
#ifdef INT64_MAX
template<>
struct FixedDataType<64, false> {
  using type = int64_t;
};
#endif
#ifdef INT128_MAX
template<>
struct FixedDataType<128, false> {
  using type = __int128_t;
};
#endif


template<size_t N, size_t M, bool fast>
struct Fixed {
public:
  static constexpr size_t bits = N;
  static constexpr size_t frac = M;
  static constexpr bool is_fast = fast;
  using intN_t = typename FixedDataType<N, fast>::type;
  using int2N_t = typename FixedDataType<2 * N, fast>::type;

  static_assert(M < N, "M must be less than N");

  static_assert(!std::is_same<intN_t, double>::value, "N is too large");
  constexpr Fixed(int v): v(v << M) {}
  constexpr Fixed(float f): v(f * (1 << M)) {}
  constexpr Fixed(double f): v(f * (1 << M)) {}
  constexpr Fixed(): v(0) {}

  static constexpr Fixed from_raw(intN_t x) {
    Fixed ret;
    ret.v = x;
    return ret;
  }

  intN_t v;

  auto operator<=>(const Fixed& other) const {
    return v <=> other.v;
  }
  bool operator==(const Fixed& other) const {
    return v == other.v;
  }

  template<size_t N2, size_t M2, bool fast2>
  constexpr operator Fixed<N2, M2, fast2>() const {
    using dst_storage = typename FixedDataType<N2, fast2>::type;

    if constexpr (M2 > M) {
      intN_t shifted = v;
      shifted <<= (M2 - M);
      return Fixed<N2, M2, fast2>::from_raw(static_cast<dst_storage>(shifted));
    } else {
      intN_t shifted = v;
      shifted >>= (M - M2);
      return Fixed<N2, M2, fast2>::from_raw(static_cast<dst_storage>(shifted));
    }
  }
  constexpr operator double() const {
    return double(v) / (1 << M);
  }
  constexpr operator float() const {
    return float(v) / (1 << M);
  }
  //  static constexpr Fixed inf = Fixed::from_raw(std::numeric_limits<intN_t>::max());
  //  static constexpr Fixed eps = Fixed::from_raw(deltas.size());
};

template<size_t N, size_t M, bool fast>
Fixed<N, M, fast> operator+(Fixed<N, M, fast> a, Fixed<N, M, fast> b) {
  return Fixed<N, M, fast>::from_raw(a.v + b.v);
}

template<size_t N, size_t M, bool fast>
Fixed<N, M, fast> operator-(Fixed<N, M, fast> a, Fixed<N, M, fast> b) {
  return Fixed<N, M, fast>::from_raw(a.v - b.v);
}

template<size_t N, size_t M, bool fast>
Fixed<N, M, fast> operator*(Fixed<N, M, fast> a, Fixed<N, M, fast> b) {
  return Fixed<N, M, fast>::from_raw((static_cast<typename Fixed<N, M, fast>::int2N_t>(a.v) * b.v) / (1 << M));
}

template<size_t N, size_t M, bool fast>
Fixed<N, M, fast> operator/(Fixed<N, M, fast> a, Fixed<N, M, fast> b) {
//  assert(b.v != 0);
  return Fixed<N, M, fast>::from_raw(((static_cast<typename Fixed<N, M, fast>::int2N_t>(a.v) * (1 << M))) / b.v);
}

template<size_t N, size_t M, bool fast>
Fixed<N, M, fast> &operator+=(Fixed<N, M, fast> &a, Fixed<N, M, fast> b) {
  return a = a + b;
}

template<size_t N, size_t M, bool fast>
Fixed<N, M, fast> &operator-=(Fixed<N, M, fast> &a, Fixed<N, M, fast> b) {
  return a = a - b;
}

template<size_t N, size_t M, bool fast>
Fixed<N, M, fast> &operator*=(Fixed<N, M, fast> &a, Fixed<N, M, fast> b) {
  return a = a * b;
}

template<size_t N, size_t M, bool fast>
Fixed<N, M, fast> &operator/=(Fixed<N, M, fast> &a, Fixed<N, M, fast> b) {
  return a = a / b;
}

template<size_t N, size_t M, bool fast>
Fixed<N, M, fast> operator-(Fixed<N, M, fast> x) {
  return Fixed<N, M, fast>::from_raw(-x.v);
}

template<size_t N, size_t M, bool fast>
Fixed<N, M, fast> abs(Fixed<N, M, fast> x) {
  if (x.v < 0) {
    x.v = -x.v;
  }
  return x;
}

template<size_t N, size_t M, bool fast>
ostream &operator<<(ostream &out, Fixed<N, M, fast> x) {
  return out << x.v / (double) (1 << M);
}

class meowfluid_base {
public:
  virtual void run() = 0;
  virtual void init(const InitArgs &initArgs) = 0;
  virtual void save() const = 0;
  virtual void load() = 0;
  virtual void setStop() = 0;
};

template<typename Tp, typename Tv, typename Tf, size_t N, size_t M>
class meowfluid : public meowfluid_base{
private:
  template<typename Fixed>
  struct VectorField {
    std::unique_ptr<std::array<std::array<array<Fixed, deltas.size()>, M>, N>> ptr_v;
    Fixed &add(int x, int y, int dx, int dy, Fixed dv) {
      return get(x, y, dx, dy) += dv;
    }

    void setZeros() {
      auto &v = *ptr_v.get();
      for (int x = 0; x < v.size(); ++x) {
        for (int y = 0; y < v[0].size(); ++y) {
          for (int d = 0; d < v[0][0].size(); ++d) {
            v[x][y][d] = Fixed(0);
          }
        }
      }
    }

    Fixed &get(int x, int y, int dx, int dy) {
      //          size_t i = ranges::find(deltas, pair(dx, dy)) - deltas.begin();
      size_t i = find(deltas.begin(), deltas.end(), pair(dx, dy)) - deltas.begin();
      assert(i < deltas.size());
      return (*ptr_v.get())[x][y][i];
    }
  };

  std::unique_ptr<std::array<std::array<int, M>, N>> ptr_dirs;
  mt19937 rnd;
  bool stop = false;
  double rho[256];
  std::unique_ptr<std::array<std::array<char, M>, N>> ptr_field;
  VectorField<Tv> velocity;
  VectorField<Tf> velocity_flow;
  std::unique_ptr<std::array<std::array<Tp, M>, N>> ptr_p, ptr_old_p;
  std::unique_ptr<std::array<std::array<int, M>, N>> ptr_last_use;
  int UT = 0;

  size_t iter = 0;
  static constexpr size_t T = 1'000;
  double g;

public:
  meowfluid() = default;
  void setStop() override { stop = true; }
  void init(const InitArgs &initArgs) override {
    ptr_dirs = std::make_unique<std::array<std::array<int,M>,N>>();
    ptr_field = std::make_unique<std::array<std::array<char,M>,N>>();
    ptr_p = std::make_unique<std::array<std::array<Tp,M>,N>>();
    ptr_old_p = std::make_unique<std::array<std::array<Tp,M>,N>>();
    ptr_last_use = std::make_unique<std::array<std::array<int,M>,N>>();
    velocity.ptr_v = std::make_unique<std::array<std::array<array<Tv, deltas.size()>,M>,N>>();
    velocity_flow.ptr_v = std::make_unique<std::array<std::array<array<Tf, deltas.size()>,M>,N>>();
    auto &field = *ptr_field.get();
    auto &dirs = *ptr_dirs.get();

    auto &tmp_field = initArgs.tmp_field;
    if (tmp_field.size() != N || tmp_field[0].size() != M) {
      std::cout << "Something went wrong, field size does not match template size" << std::endl;
      exit(1);
      return;
    }
    iter=0;
    rnd = mt19937(1337);
    g = initArgs.G;
    rho[int(' ')] = initArgs.RHO_A;
    rho[int('.')] = initArgs.RHO_F;
    for (size_t x = 0; x < N; ++x) {
      for (size_t y = 0; y < M; ++y) {
        field[x][y] = tmp_field[x][y];
        if (field[x][y] == '#')
          continue;
        for (auto [dx, dy] : deltas) {
          dirs[x][y] += (field[x + dx][y + dy] != '#');
        }
      }
    }
  }

  ~meowfluid() = default;

  void run() override {
    auto &field = *ptr_field.get();
    auto &dirs = *ptr_dirs.get();
    auto &p = *ptr_p.get();
    auto &old_p = *ptr_old_p.get();
    auto &last_use = *ptr_last_use.get();

    for (; iter < T && !stop; ++iter) {
//      Tp total_delta_p = 0;
      // Apply external forces
      for (size_t x = 0; x < N; ++x) {
        for (size_t y = 0; y < M; ++y) {
          if (field[x][y] == '#')
            continue;
          if (field[x + 1][y] != '#')
            velocity.add(x, y, 1, 0, g);
        }
      }

      // Apply forces from p
      old_p = p;
      //      memcpy(old_p, p, sizeof(p));
      for (size_t x = 0; x < N; ++x) {
        for (size_t y = 0; y < M; ++y) {
          if (field[x][y] == '#')
            continue;
          for (auto [dx, dy] : deltas) {
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] != '#' && old_p[nx][ny] < old_p[x][y]) {
              auto delta_p = old_p[x][y] - old_p[nx][ny];
              auto force = delta_p;
              auto &contr = velocity.get(nx, ny, -dx, -dy);
              if (contr * static_cast<Tv>(rho[(int) field[nx][ny]]) >= static_cast<Tv>(force)) {
                contr -= static_cast<Tv>(force) / static_cast<Tv>(rho[(int) field[nx][ny]]);
                continue;
              }
              force -= static_cast<Tp>(contr) * static_cast<Tp>(rho[(int) field[nx][ny]]);
              contr = 0;
              velocity.add(x, y, dx, dy, force / static_cast<Tp>(rho[(int) field[x][y]]));
              p[x][y] -= force / Tp(dirs[x][y]);
//              total_delta_p -= force / Fixed(dirs[x][y]);
            }
          }
        }
      }

      // Make flow from velocities
//      velocity_flow = {};
      velocity_flow.setZeros();
      bool prop = false;
      do {
        UT += 2;
        prop = 0;
        for (size_t x = 0; x < N; ++x) {
          for (size_t y = 0; y < M; ++y) {
            if (field[x][y] != '#' && last_use[x][y] != UT) {
              auto [t, local_prop, _] = propagate_flow(x, y, 1);
              if (t > static_cast<Tv>(0)) {
                prop = 1;
              }
            }
          }
        }
      } while (prop);

      // Recalculate p with kinetic energy
      for (size_t x = 0; x < N; ++x) {
        for (size_t y = 0; y < M; ++y) {
          if (field[x][y] == '#')
            continue;
          for (auto [dx, dy] : deltas) {
            auto old_v = velocity.get(x, y, dx, dy);
            auto new_v = velocity_flow.get(x, y, dx, dy);
            if (old_v > static_cast<Tv>(0)) {
              assert(static_cast<Tv>(new_v) <= old_v);
              velocity.get(x, y, dx, dy) = static_cast<Tv>(new_v);
              auto force = static_cast<Tp>((old_v - static_cast<Tv>(new_v))) * static_cast<Tp>(rho[(int) field[x][y]]);
              if (field[x][y] == '.')
                force *= Tp(0.8);
              if (field[x + dx][y + dy] == '#') {
                p[x][y] += force / Tp(dirs[x][y]);
//                total_delta_p += force / Fixed(dirs[x][y]);
              } else {
                p[x + dx][y + dy] += force / Tp(dirs[x + dx][y + dy]);
//                total_delta_p += force / Fixed(dirs[x + dx][y + dy]);
              }
            }
          }
        }
      }

      UT += 2;
      prop = false;
      for (size_t x = 0; x < N; ++x) {
        for (size_t y = 0; y < M; ++y) {
          if (field[x][y] != '#' && last_use[x][y] != UT) {
            if (random01<Tv>() < move_prob(x, y)) {
              prop = true;
              propagate_move(x, y, true);
            } else {
              propagate_stop(x, y, true);
            }
          }
        }
      }

      if (prop || iter % 1000 == 0) {
        system("clear");
        cout << "Tick " << iter << ":\n";
        for (size_t x = 0; x < field.size(); ++x) {
          for (size_t y = 0; y < field[0].size(); ++y) {
            cout << field[x][y];
          }
          cout << "\n";
        }
        std::cout << std::endl;
      }
    }
    save();
  }
  struct ParticleParams {
    char type;
    Tp cur_p;
    array<Tv, deltas.size()> v;
  };

  void swap_with(int x, int y, ParticleParams &pp) {

    auto &field = *ptr_field.get();
    auto &p = *ptr_p.get();
    auto &v = *(velocity.ptr_v.get());
    swap(field[x][y], pp.type);
    swap(p[x][y], pp.cur_p);
    swap(v[x][y], pp.v);
  }

  bool propagate_move(int x, int y, bool is_first) {
    auto &field = *ptr_field.get();
    auto &last_use = *ptr_last_use.get();
    last_use[x][y] = UT - is_first;
    bool ret = false;
    int nx = -1, ny = -1;
    do {
      std::array<Tv, deltas.size()> tres;
      Tv sum = 0;
      for (size_t i = 0; i < deltas.size(); ++i) {
        auto [dx, dy] = deltas[i];
        int nx = x + dx, ny = y + dy;
        if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
          tres[i] = sum;
          continue;
        }
        auto v = velocity.get(x, y, dx, dy);
        if (v < static_cast<Tv>(0)) {
          tres[i] = sum;
          continue;
        }
        sum += v;
        tres[i] = sum;
      }

      if (sum == static_cast<Tv>(0)) {
        break;
      }

      Tp p = random01<Tp>() * static_cast<Tp>(sum);
      size_t d = std::upper_bound(tres.begin(), tres.end(), static_cast<Tv>(p)) - tres.begin();

      auto [dx, dy] = deltas[d];
      nx = x + dx;
      ny = y + dy;
      assert(velocity.get(x, y, dx, dy) > static_cast<Tv>(0) && field[nx][ny] != '#' && last_use[nx][ny] < UT);

      ret = (last_use[nx][ny] == UT - 1 || propagate_move(nx, ny, false));
    } while (!ret);
    last_use[x][y] = UT;
    for (size_t i = 0; i < deltas.size(); ++i) {
      auto [dx, dy] = deltas[i];
      int nx = x + dx, ny = y + dy;
      if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) < static_cast<Tv>(0)) {
        propagate_stop(nx, ny);
      }
    }
    if (ret) {
      if (!is_first) {
        ParticleParams pp{};
        swap_with(x, y, pp);
        swap_with(nx, ny, pp);
        swap_with(x, y, pp);
      }
    }
    return ret;
  }
  Tv move_prob(int x, int y) {
    auto &field = *ptr_field.get();
    auto &last_use = *ptr_last_use.get();
    Tv sum = 0;
    for (size_t i = 0; i < deltas.size(); ++i) {
      auto [dx, dy] = deltas[i];
      int nx = x + dx, ny = y + dy;
      if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
        continue;
      }
      auto v = velocity.get(x, y, dx, dy);
      if (v < static_cast<Tv>(0)) {
        continue;
      }
      sum += v;
    }
    return sum;
  }
  void propagate_stop(int x, int y, bool force = false) {
    auto &field = *ptr_field.get();
    auto &last_use = *ptr_last_use.get();
    if (!force) {
      bool stop = true;
      for (auto [dx, dy] : deltas) {
        int nx = x + dx, ny = y + dy;
        if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) > static_cast<Tv>(0)) {
          stop = false;
          break;
        }
      }
      if (!stop) {
        return;
      }
    }
    last_use[x][y] = UT;
    for (auto [dx, dy] : deltas) {
      int nx = x + dx, ny = y + dy;
      if (field[nx][ny] == '#' || last_use[nx][ny] == UT || velocity.get(x, y, dx, dy) > static_cast<Tv>(0)) {
        continue;
      }
      propagate_stop(nx, ny);
    }
  }
  template<typename Fixed>
  Fixed random01() {
    if constexpr (std::is_same_v<Fixed, double>) {
      return (rnd() / (double)rnd.max());
    }
    else if constexpr (std::is_same_v<Fixed, float>) {
      return (rnd() / (float)rnd.max());
    }
    else {
      Fixed tmp = Fixed::from_raw((rnd() & ((1 << Fixed::frac) - 1)));
      assert(tmp >= static_cast<Fixed>(0) && tmp <= static_cast<Fixed>(1));
      return tmp;
    }
  }
  tuple<Tv, bool, pair<int, int>> propagate_flow(int x, int y, Tv lim) {
    auto &field = *ptr_field.get();
    auto &last_use = *ptr_last_use.get();
    last_use[x][y] = UT - 1;
    Tv ret = 0;
    for (auto [dx, dy] : deltas) {
      int nx = x + dx, ny = y + dy;
      if (field[nx][ny] != '#' && last_use[nx][ny] < UT) {
        auto cap = velocity.get(x, y, dx, dy);
        auto flow = velocity_flow.get(x, y, dx, dy);
        if (flow == static_cast<Tf>(cap)) {
          continue;
        }
        // assert(v >= velocity_flow.get(x, y, dx, dy));
        auto vp = min(lim, cap - static_cast<Tv>(flow));
        if (last_use[nx][ny] == UT - 1) {
          velocity_flow.add(x, y, dx, dy, vp);
          last_use[x][y] = UT;
          // cerr << x << " " << y << " -> " << nx << " " << ny << " " << vp << " / " << lim << "\n";
          return {vp, 1, {nx, ny}};
        }
        auto [t, prop, end] = propagate_flow(nx, ny, vp);
        ret += t;
        if (prop) {
          velocity_flow.add(x, y, dx, dy, t);
          last_use[x][y] = UT;
          // cerr << x << " " << y << " -> " << nx << " " << ny << " " << t << " / " << lim << "\n";
          return {t, prop && end != pair(x, y), end};
        }
      }
    }
    last_use[x][y] = UT;
    return {ret, 0, {0, 0}};
  }

  void save() const override {

    ofstream out;
    out.open("meowfluid.save", std::ios::binary | std::ios::out);

    auto &field = *ptr_field.get();
    auto &dirs = *ptr_dirs.get();
    auto &p = *ptr_p.get();
    auto &old_p = *ptr_old_p.get();
    auto &last_use = *ptr_last_use.get();
    auto &v = *velocity.ptr_v.get();
    auto &v_flow = *velocity_flow.ptr_v.get();
    size_t NN = N;
    size_t MM = M;
    out.write(reinterpret_cast<const char*>(&NN), sizeof(NN));
    out.write(reinterpret_cast<const char*>(&MM), sizeof(MM));
    out.write(reinterpret_cast<const char*>(&iter), sizeof(size_t));
    out.write(reinterpret_cast<const char*>(&UT), sizeof(UT));
    out.write(reinterpret_cast<const char*>(&g), sizeof(g));


    out.write(reinterpret_cast<const char*>(&rho[int('.')]),sizeof(double));
    out.write(reinterpret_cast<const char*>(&rho[int(' ')]),sizeof(double));

    for (size_t i = 0; i < N; ++i) {
      out.write(reinterpret_cast<const char*>(field[i].data()), M * sizeof(char));
      out.write(reinterpret_cast<const char*>(dirs[i].data()), M * sizeof(int));
      out.write(reinterpret_cast<const char*>(p[i].data()), M * sizeof(Tp));

      out.write(reinterpret_cast<const char*>(old_p[i].data()), M * sizeof(Tp));
      out.write(reinterpret_cast<const char*>(last_use[i].data()), M * sizeof(int));
    }

    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < M; ++j) {
        out.write(reinterpret_cast<const char*>(v[i][j].data()), deltas.size() * sizeof(Tv));
        out.write(reinterpret_cast<const char*>(v_flow[i][j].data()), deltas.size() * sizeof(Tf));
      }
    }

    out.close();
  }

  void load() override {
    fstream in("meowfluid.save", std::ios::binary | std::ios::in);

    ptr_dirs = std::make_unique<std::array<std::array<int,M>,N>>();
    ptr_field = std::make_unique<std::array<std::array<char,M>,N>>();
    ptr_p = std::make_unique<std::array<std::array<Tp,M>,N>>();
    ptr_old_p = std::make_unique<std::array<std::array<Tp,M>,N>>();
    ptr_last_use = std::make_unique<std::array<std::array<int,M>,N>>();
    velocity.ptr_v = std::make_unique<std::array<std::array<array<Tv, deltas.size()>,M>,N>>();
    velocity_flow.ptr_v = std::make_unique<std::array<std::array<array<Tf, deltas.size()>,M>,N>>();

    auto &field = *ptr_field.get();
    auto &dirs = *ptr_dirs.get();
    auto &p = *ptr_p.get();
    auto &old_p = *ptr_old_p.get();
    auto &last_use = *ptr_last_use.get();
    auto &v = *velocity.ptr_v.get();
    auto &v_flow = *velocity_flow.ptr_v.get();

    size_t NN, MM;
    in.read(reinterpret_cast<char*>(&NN), sizeof(NN));
    in.read(reinterpret_cast<char*>(&MM), sizeof(MM));

    if (NN != N || MM != M) {
      std::cout << "Bad size for load" << std::endl;
      exit(0);
      return;
    }

    in.read(reinterpret_cast<char*>(&iter), sizeof(iter));
    in.read(reinterpret_cast<char*>(&UT), sizeof(UT));
    in.read(reinterpret_cast<char*>(&g), sizeof(g));


    in.read(reinterpret_cast<char*>(&rho[int('.')]),sizeof(double));
    in.read(reinterpret_cast<char*>(&rho[int(' ')]),sizeof(double));

    for (size_t i = 0; i < N; ++i) {
      in.read(reinterpret_cast<char*>(field[i].data()), M * sizeof(char));
      in.read(reinterpret_cast<char*>(dirs[i].data()), M * sizeof(int));
      in.read(reinterpret_cast<char*>(p[i].data()), M * sizeof(Tp));

      in.read(reinterpret_cast<char*>(old_p[i].data()), M * sizeof(Tp));
      in.read(reinterpret_cast<char*>(last_use[i].data()), M * sizeof(int));
    }

    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < M; ++j) {
        in.read(reinterpret_cast<char*>(v[i][j].data()), deltas.size() * sizeof(Tv));
        in.read(reinterpret_cast<char*>(v_flow[i][j].data()), deltas.size() * sizeof(Tf));
      }
    }

    in.close();

    rnd = mt19937(1337);
  }
};


template<typename Tp, typename Tv, typename Tf>
class meowfluidDynamic : public meowfluid_base{
private:

  template<typename Fixed>
  struct VectorField {
    std::vector<std::vector<array<Fixed, deltas.size()>>> v;
    Fixed &add(int x, int y, int dx, int dy, Fixed dv) {
      return get(x, y, dx, dy) += dv;
    }

    void setZeros() {
      for (int x = 0; x < v.size(); ++x) {
        for (int y = 0; y < v[0].size(); ++y) {
          for (int d = 0; d < v[0][0].size(); ++d) {
            v[x][y][d] = Fixed(0);
          }
        }
      }
    }

    Fixed &get(int x, int y, int dx, int dy) {
      //          size_t i = ranges::find(deltas, pair(dx, dy)) - deltas.begin();
      size_t i = find(deltas.begin(), deltas.end(), pair(dx, dy)) - deltas.begin();
      assert(i < deltas.size());
      return v[x][y][i];
    }
  };

  std::vector<std::vector<int>> dirs;
  mt19937 rnd;
  bool stop = false;
  double rho[256];
  std::vector<std::vector<char>> field;
  VectorField<Tv> velocity;
  VectorField<Tf> velocity_flow;
  std::vector<std::vector<Tp>> p, old_p;
  std::vector<std::vector<int>> last_use;
  int UT = 0;
  static constexpr size_t T = 1'000;
  double g = 0.1;
  size_t N, M;
  size_t iter = 0;

public:
  meowfluidDynamic() = default;
  void setStop() override { stop = true; }
  void init(const InitArgs &initArgs) override {
    auto &tmp_field = initArgs.tmp_field;
    N = tmp_field.size(), M = tmp_field[0].size();
    dirs = std::vector<std::vector<int>>(tmp_field.size(), std::vector<int>(tmp_field[0].size()));
    field = std::vector<std::vector<char>>(tmp_field.size(), std::vector<char>(tmp_field[0].size()));
    p = std::vector<std::vector<Tp>>(tmp_field.size(), std::vector<Tp>(tmp_field[0].size()));
    old_p = std::vector<std::vector<Tp>>(tmp_field.size(), std::vector<Tp>(tmp_field[0].size()));
    last_use = std::vector<std::vector<int>>(tmp_field.size(), std::vector<int>(tmp_field[0].size()));
    velocity.v = std::vector<std::vector<array<Tv, deltas.size()>>>(tmp_field.size(), std::vector<array<Tv, deltas.size()>>(tmp_field[0].size(), array<Tv, deltas.size()>()));
    velocity_flow.v = std::vector<std::vector<array<Tf, deltas.size()>>>(tmp_field.size(), std::vector<array<Tf, deltas.size()>>(tmp_field[0].size(), array<Tf, deltas.size()>()));

    rnd = mt19937(1337);
    iter=0;
    g = initArgs.G;
    rho[int(' ')] = initArgs.RHO_A;
    rho[int('.')] = initArgs.RHO_F;
    for (size_t x = 0; x < N; ++x) {
      for (size_t y = 0; y < M; ++y) {
        field[x][y] = tmp_field[x][y];
        if (field[x][y] == '#')
          continue;
        for (auto [dx, dy] : deltas) {
          dirs[x][y] += (field[x + dx][y + dy] != '#');
        }
      }
    }
  }

  ~meowfluidDynamic() = default;

  void run() override {

    for (; iter < T && !stop; ++iter) {
//      Tp total_delta_p = 0;
      // Apply external forces
      for (size_t x = 0; x < N; ++x) {
        for (size_t y = 0; y < M; ++y) {
          if (field[x][y] == '#')
            continue;
          if (field[x + 1][y] != '#')
            velocity.add(x, y, 1, 0, g);
        }
      }

      // Apply forces from p
      old_p = p;
      //      memcpy(old_p, p, sizeof(p));
      for (size_t x = 0; x < N; ++x) {
        for (size_t y = 0; y < M; ++y) {
          if (field[x][y] == '#')
            continue;
          for (auto [dx, dy] : deltas) {
            int nx = x + dx, ny = y + dy;
            if (field[nx][ny] != '#' && old_p[nx][ny] < old_p[x][y]) {
              auto delta_p = old_p[x][y] - old_p[nx][ny];
              auto force = delta_p;
              auto &contr = velocity.get(nx, ny, -dx, -dy);
              if (contr * static_cast<Tv>(rho[(int) field[nx][ny]]) >= static_cast<Tv>(force)) {
                contr -= static_cast<Tv>(force) / static_cast<Tv>(rho[(int) field[nx][ny]]);
                continue;
              }
              force -= static_cast<Tp>(contr) * static_cast<Tp>(rho[(int) field[nx][ny]]);
              contr = 0;
              velocity.add(x, y, dx, dy, force / static_cast<Tp>(rho[(int) field[x][y]]));
              p[x][y] -= force / Tp(dirs[x][y]);
//              total_delta_p -= force / Fixed(dirs[x][y]);
            }
          }
        }
      }

      // Make flow from velocities
//      velocity_flow = {};
      velocity_flow.setZeros();
      bool prop = false;
      do {
        UT += 2;
        prop = 0;
        for (size_t x = 0; x < N; ++x) {
          for (size_t y = 0; y < M; ++y) {
            if (field[x][y] != '#' && last_use[x][y] != UT) {
              auto [t, local_prop, _] = propagate_flow(x, y, 1);
              if (t > static_cast<Tv>(0)) {
                prop = 1;
              }
            }
          }
        }
      } while (prop);

      // Recalculate p with kinetic energy
      for (size_t x = 0; x < N; ++x) {
        for (size_t y = 0; y < M; ++y) {
          if (field[x][y] == '#')
            continue;
          for (auto [dx, dy] : deltas) {
            auto old_v = velocity.get(x, y, dx, dy);
            auto new_v = velocity_flow.get(x, y, dx, dy);
            if (old_v > static_cast<Tv>(0)) {
              assert(static_cast<Tv>(new_v) <= old_v);
              velocity.get(x, y, dx, dy) = static_cast<Tv>(new_v);
              auto force = static_cast<Tp>((old_v - static_cast<Tv>(new_v))) * static_cast<Tp>(rho[(int) field[x][y]]);
              if (field[x][y] == '.')
                force *= Tp(0.8);
              if (field[x + dx][y + dy] == '#') {
                p[x][y] += force / Tp(dirs[x][y]);
//                total_delta_p += force / Fixed(dirs[x][y]);
              } else {
                p[x + dx][y + dy] += force / Tp(dirs[x + dx][y + dy]);
//                total_delta_p += force / Fixed(dirs[x + dx][y + dy]);
              }
            }
          }
        }
      }

      UT += 2;
      prop = false;
      for (size_t x = 0; x < N; ++x) {
        for (size_t y = 0; y < M; ++y) {
          if (field[x][y] != '#' && last_use[x][y] != UT) {
            if (random01<Tv>() < move_prob(x, y)) {
              prop = true;
              propagate_move(x, y, true);
            } else {
              propagate_stop(x, y, true);
            }
          }
        }
      }

      if (prop || iter % 1000 == 0) {
        system("clear");
        cout << "Tick " << iter << ":\n";
        for (size_t x = 0; x < field.size(); ++x) {
          for (size_t y = 0; y < field[0].size(); ++y) {
            cout << field[x][y];
          }
          cout << "\n";
        }
        std::cout << std::endl;
      }
    }
    save();
  }
  struct ParticleParams {
    char type;
    Tp cur_p;
    array<Tv, deltas.size()> v;
  };

  void swap_with(int x, int y, ParticleParams &pp) {
    swap(field[x][y], pp.type);
    swap(p[x][y], pp.cur_p);
    swap(velocity.v[x][y], pp.v);
  }

  bool propagate_move(int x, int y, bool is_first) {
    last_use[x][y] = UT - is_first;
    bool ret = false;
    int nx = -1, ny = -1;
    do {
      std::array<Tv, deltas.size()> tres;
      Tv sum = 0;
      for (size_t i = 0; i < deltas.size(); ++i) {
        auto [dx, dy] = deltas[i];
        int nx = x + dx, ny = y + dy;
        if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
          tres[i] = sum;
          continue;
        }
        auto v = velocity.get(x, y, dx, dy);
        if (v < static_cast<Tv>(0)) {
          tres[i] = sum;
          continue;
        }
        sum += v;
        tres[i] = sum;
      }

      if (sum == static_cast<Tv>(0)) {
        break;
      }

      Tv p = random01<Tv>() * sum;
      size_t d = std::upper_bound(tres.begin(), tres.end(), p) - tres.begin();

      auto [dx, dy] = deltas[d];
      nx = x + dx;
      ny = y + dy;
      assert(velocity.get(x, y, dx, dy) > static_cast<Tv>(0) && field[nx][ny] != '#' && last_use[nx][ny] < UT);

      ret = (last_use[nx][ny] == UT - 1 || propagate_move(nx, ny, false));
    } while (!ret);
    last_use[x][y] = UT;
    for (size_t i = 0; i < deltas.size(); ++i) {
      auto [dx, dy] = deltas[i];
      int nx = x + dx, ny = y + dy;
      if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) < static_cast<Tv>(0)) {
        propagate_stop(nx, ny);
      }
    }
    if (ret) {
      if (!is_first) {
        ParticleParams pp{};
        swap_with(x, y, pp);
        swap_with(nx, ny, pp);
        swap_with(x, y, pp);
      }
    }
    return ret;
  }

  Tv move_prob(int x, int y) {
    Tv sum = 0;
    for (size_t i = 0; i < deltas.size(); ++i) {
      auto [dx, dy] = deltas[i];
      int nx = x + dx, ny = y + dy;
      if (field[nx][ny] == '#' || last_use[nx][ny] == UT) {
        continue;
      }
      auto v = velocity.get(x, y, dx, dy);
      if (v < static_cast<Tv>(0)) {
        continue;
      }
      sum += v;
    }
    return sum;
  }

  void propagate_stop(int x, int y, bool force = false) {
    if (!force) {
      bool stop = true;
      for (auto [dx, dy] : deltas) {
        int nx = x + dx, ny = y + dy;
        if (field[nx][ny] != '#' && last_use[nx][ny] < UT - 1 && velocity.get(x, y, dx, dy) > static_cast<Tv>(0)) {
          stop = false;
          break;
        }
      }
      if (!stop) {
        return;
      }
    }
    last_use[x][y] = UT;
    for (auto [dx, dy] : deltas) {
      int nx = x + dx, ny = y + dy;
      if (field[nx][ny] == '#' || last_use[nx][ny] == UT || velocity.get(x, y, dx, dy) > static_cast<Tv>(0)) {
        continue;
      }
      propagate_stop(nx, ny);
    }
  }
  template<typename Fixed>
  Fixed random01() {
    if constexpr (std::is_same_v<Fixed, double>) {
      return (rnd() / (double)rnd.max());
    }
    else if constexpr (std::is_same_v<Fixed, float>) {
      return (rnd() / (float)rnd.max());
    }
    else {
      Fixed tmp = Fixed::from_raw((rnd() & ((1 << Fixed::frac) - 1)));
      assert(tmp >= static_cast<Fixed>(0) && tmp <= static_cast<Fixed>(1));
      return tmp;
    }
  }
  tuple<Tv, bool, pair<int, int>> propagate_flow(int x, int y, Tv lim) {
    last_use[x][y] = UT - 1;
    Tv ret = 0;
    for (auto [dx, dy] : deltas) {
      int nx = x + dx, ny = y + dy;
      if (field[nx][ny] != '#' && last_use[nx][ny] < UT) {
        auto cap = velocity.get(x, y, dx, dy);
        auto flow = velocity_flow.get(x, y, dx, dy);
        if (flow == static_cast<Tf>(cap)) {
          continue;
        }
        // assert(v >= velocity_flow.get(x, y, dx, dy));
        auto vp = min(lim, cap - static_cast<Tv>(flow));
        if (last_use[nx][ny] == UT - 1) {
          velocity_flow.add(x, y, dx, dy, vp);
          last_use[x][y] = UT;
          // cerr << x << " " << y << " -> " << nx << " " << ny << " " << vp << " / " << lim << "\n";
          return {vp, 1, {nx, ny}};
        }
        auto [t, prop, end] = propagate_flow(nx, ny, vp);
        ret += t;
        if (prop) {
          velocity_flow.add(x, y, dx, dy, t);
          last_use[x][y] = UT;
          // cerr << x << " " << y << " -> " << nx << " " << ny << " " << t << " / " << lim << "\n";
          return {t, prop && end != pair(x, y), end};
        }
      }
    }
    last_use[x][y] = UT;
    return {ret, 0, {0, 0}};
  }

  void save() const override {

    ofstream out;
    out.open("meowfluid.save", std::ios::binary | std::ios::out);

    out.write(reinterpret_cast<const char*>(&N), sizeof(N));
    out.write(reinterpret_cast<const char*>(&M), sizeof(M));
    out.write(reinterpret_cast<const char*>(&iter), sizeof(iter));
    out.write(reinterpret_cast<const char*>(&UT), sizeof(UT));
    out.write(reinterpret_cast<const char*>(&g), sizeof(g));


    out.write(reinterpret_cast<const char*>(&rho[int('.')]),sizeof(double));
    out.write(reinterpret_cast<const char*>(&rho[int(' ')]),sizeof(double));

    for (size_t i = 0; i < N; ++i) {
      out.write(reinterpret_cast<const char*>(field[i].data()), M * sizeof(char));
      out.write(reinterpret_cast<const char*>(dirs[i].data()), M * sizeof(int));
      out.write(reinterpret_cast<const char*>(p[i].data()), M * sizeof(Tp));

      out.write(reinterpret_cast<const char*>(old_p[i].data()), M * sizeof(Tp));
      out.write(reinterpret_cast<const char*>(last_use[i].data()), M * sizeof(int));
    }

    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < M; ++j) {
        out.write(reinterpret_cast<const char*>(velocity.v[i][j].data()), deltas.size() * sizeof(Tv));
        out.write(reinterpret_cast<const char*>(velocity_flow.v[i][j].data()), deltas.size() * sizeof(Tf));
      }
    }

    out.close();
  }

  void load() override {
    fstream in("meowfluid.save", std::ios::binary | std::ios::in);

    in.read(reinterpret_cast<char*>(&N), sizeof(N));
    in.read(reinterpret_cast<char*>(&M), sizeof(M));

    dirs = std::vector<std::vector<int>>(N, std::vector<int>(M));
    field = std::vector<std::vector<char>>(N, std::vector<char>(M));
    p = std::vector<std::vector<Tp>>(N, std::vector<Tp>(M));
    old_p = std::vector<std::vector<Tp>>(N, std::vector<Tp>(M));
    last_use = std::vector<std::vector<int>>(N, std::vector<int>(M));
    velocity.v = std::vector<std::vector<array<Tv, deltas.size()>>>(N, std::vector<array<Tv, deltas.size()>>(M, array<Tv, deltas.size()>()));
    velocity_flow.v = std::vector<std::vector<array<Tf, deltas.size()>>>(N, std::vector<array<Tf, deltas.size()>>(M, array<Tf, deltas.size()>()));

    rnd = mt19937(1337);

    in.read(reinterpret_cast<char*>(&iter), sizeof(iter));
    in.read(reinterpret_cast<char*>(&UT), sizeof(UT));
    in.read(reinterpret_cast<char*>(&g), sizeof(g));


    in.read(reinterpret_cast<char*>(&rho[int('.')]),sizeof(double));
    in.read(reinterpret_cast<char*>(&rho[int(' ')]),sizeof(double));

    for (size_t i = 0; i < N; ++i) {
      in.read(reinterpret_cast<char*>(field[i].data()), M * sizeof(char));
      in.read(reinterpret_cast<char*>(dirs[i].data()), M * sizeof(int));
      in.read(reinterpret_cast<char*>(p[i].data()), M * sizeof(Tp));

      in.read(reinterpret_cast<char*>(old_p[i].data()), M * sizeof(Tp));
      in.read(reinterpret_cast<char*>(last_use[i].data()), M * sizeof(int));
    }

    for (size_t i = 0; i < N; ++i) {
      for (size_t j = 0; j < M; ++j) {
        in.read(reinterpret_cast<char*>(velocity.v[i][j].data()), deltas.size() * sizeof(Tv));
        in.read(reinterpret_cast<char*>(velocity_flow.v[i][j].data()), deltas.size() * sizeof(Tf));
      }
    }

    in.close();
  }
};


static std::unordered_map<std::string, std::function<void(int,int)>> g_registry;

static std::vector<std::string> g_registry_types;

static meowfluid_base *mf = nullptr;

template <typename T>
std::string getTypeName() {
  if constexpr (T::is_fast) {
    return "Fast_Fixed("+std::to_string(T::bits)+","+std::to_string(T::frac)+")";
  }
  return "Fixed("+std::to_string(T::bits)+","+std::to_string(T::frac)+")";
}

template <>
std::string getTypeName<double>() {
  return "double";
}
template <>
std::string getTypeName<float>() {
  return "float";
}

template <typename Tp, typename Tv, typename Tf>
std::string getTypeNames() {
  return getTypeName<Tp>() +"|"+getTypeName<Tv>() +"|"+getTypeName<Tf>();
}

// Макросы для размеров
#define S(N,M) HANDLE_SIZE(N,M)


template<typename Tp, typename Tv, typename Tf, int N, int M>
void registryMeowfluid(std::unordered_map<std::string, std::function<void(int,int)>> &g_registry, int w, int h,
                       const InitArgs &initArgs) {

  auto typeString = getTypeNames<Tp, Tv, Tf>()+"|"+std::to_string(N)+"|"+std::to_string(M);

  g_registry[typeString] = [](int w, int h) {
    if (w == N && h == M) {
      mf = new meowfluid<Tp, Tv, Tf, N, M>();
    }
  };
}


#ifdef TYPES
#define FLOAT float
#define DOUBLE double
#define FIXED(n, m) Fixed<n, m, false>
#define FAST_FIXED(n, m) Fixed<n, m, true>
#else
static_assert(false, "There is no TYPES");
#endif


template <typename Tp, typename Tv, typename Tf>
void test_7() {
//  std::cout << "test_7: " << getTypeName<Tp>() << " " << getTypeName<Tv>() << " " << getTypeName<Tf>() << std::endl;
  int tmp;
#ifdef SIZES
#define HANDLE_SIZE(N,M) registryMeowfluid<Tp, Tv, Tf, N, M>(g_registry, N, M, initArgs);tmp=5
  SIZES;
#undef HANDLE_SIZE
#endif

  auto typeString = getTypeNames<Tp, Tv, Tf>()+"|dynamic";

  g_registry[typeString] = [](int w, int h) {
    mf = new meowfluidDynamic<Tp, Tv, Tf>();
  };
}

template <typename T,typename T2, typename... _Types>
void test_6() {
  (test_7<T, T2, _Types>(), ...);
}

template <typename T, typename T2>
void test_5() {
  test_6<T, T2, TYPES>();
}

template <typename T, typename... _Types>
void test_4() {
  (test_5<T, _Types>(), ...);
}

template <typename T>
void test_3() {
  test_4<T, TYPES>();
}

//template <typename T>
//void test_2() {
//  int tmp;
//#ifdef SIZES
//#define HANDLE_SIZE(N,M) registryMeowfluid<T, N, M>(g_registry, N, M, initArgs);tmp=5
//  SIZES;
//#undef HANDLE_SIZE
//#endif
//
//  auto typeString = getTypeName<T>()+"|dynamic";
//
//  g_registry[typeString] = [](int w, int h) {
//
//    mf = new meowfluidDynamic<T>();
//  };
//}

template <typename... _Types>
void test() {
//  (test_2<_Types>(), ...);

  (test_3<_Types>(), ...);
}

extern "C" void sigint_handler(int signum) {
  if (mf) {
    mf->setStop();
  }
}

int main(int argc, char** argv) {
  test<TYPES>();

  std::string pType, vType, vFlowType, inputFile;
  bool is_load = false;

  for(int i=1; i<argc; ++i) {
    std::string arg = argv[i];
    if(arg.rfind("--p-type=", 0) == 0) {
      pType = arg.substr(9);
    } else if(arg.rfind("--v-type=", 0) == 0) {
      vType = arg.substr(9);
    } else if(arg.rfind("--v-flow-type=", 0) == 0) {
      vFlowType = arg.substr(14);
    } else if(arg.rfind("--input=",0) == 0) {
      inputFile = arg.substr(8);
    } else if(arg == "--load") {
      is_load = true;
   }
  }

  if (is_load == (inputFile.size() > 0)) {
    std::cout << "Either load or input" << std::endl;
    return 0;
  }
  size_t N, M;

  if (inputFile.size()) {
    fstream in(inputFile, std::ios::in);
  //  double G, RHO_F, RHO_A;
    in >> initArgs.G >> initArgs.RHO_A >> initArgs.RHO_F;

    in >> N >> M;
    initArgs.tmp_field = std::vector<std::vector<char>>(N, std::vector<char>(M));
    char garbage;
    in >> std::noskipws >> garbage;
    for (int i = 0; i < N; ++i) {
      for (int j = 0; j < M; ++j) {
        in >> std::noskipws >> initArgs.tmp_field[i][j];
      }
      in >> std::noskipws >> garbage;
    }
    in.close();
  } else {
    fstream in("meowfluid.save", std::ios::in);
    in.read(reinterpret_cast<char*>(&N), sizeof(N));
    in.read(reinterpret_cast<char*>(&M), sizeof(M));
    in.close();
  }


  int w = N, h = M;

  bool handled = false;

  auto key = pType+"|"+vType+"|"+vFlowType+"|"+std::to_string(w)+"|"+std::to_string(h);

  if (g_registry.find(key) != g_registry.end()) {
    // Если нашли — вызываем функцию. Если размер подходит, будет "HANDLED".
    g_registry[key](w,h);
    handled = true;
  }

  if (!handled) {
    key = pType+"|"+vType+"|"+vFlowType+"|dynamic";

    if (g_registry.find(key) != g_registry.end()) {
      g_registry[key](w,h);
    } else {
      std::cout << "Combination of types " << pType+"|"+vType+"|"+vFlowType << " is not supported" << std::endl;
      exit(0);
    }
  }

  if (!mf) {
      std::cout << "Pupupu" << std::endl;
  }

  signal(SIGINT, &sigint_handler);

  double t;
  t = clock();
  if (is_load) {
    mf->load();
  } else {
    mf->init(initArgs);
  }
  mf->run();

  std::cout << key << "time = " << clock() - t << std::endl;

  return 0;
}
