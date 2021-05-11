#ifndef L1SEEDCONEPFJETEMULATOR_H
#define L1SEEDCONEPFJETEMULATOR_H

#include "ap_int.h"
#include "ap_fixed.h"
#include "../dataformats/layer1_emulator.h"
#include "../dataformats/jets.h"
#include <iostream>
#include <vector>
#include <numeric>
#include <algorithm>

namespace l1ct{
    class L1SCEmulatedJet : public Jet {
      public:
      std::vector<PuppiObjEmu> constituents;
    };
};

class L1SCJetEmu {
  public:
  // Data types and constants used in the FPGA and FPGA-optimized functions
  // This header file is also for use in the standalone FPGA-tools simulation
  // and thus contains no CMSSW/EDM specific content
  typedef l1ct::pt_t pt_t;
  typedef l1ct::glbeta_t etaphi_t;     // Type for eta & phi
  typedef ap_int<13> detaphi_t;        // Type for deta & dphi
  typedef ap_fixed<18,23> detaphi2_t;  // Type for deta^2 & dphi^2
  typedef ap_fixed<22,22> pt_etaphi_t; // Type for product of pt with deta & dphi
  typedef l1ct::PuppiObjEmu Particle;
  typedef l1ct::L1SCEmulatedJet Jet;

  private:
  // Configuration settings
  bool _debug;
  float _coneSize;
  unsigned _nJets;
  detaphi2_t rCone2;

  // constants for the axis update
  typedef ap_ufixed<18, -2> inv_pt_t;
  static constexpr int N_table_inv_pt = 1024;
  inv_pt_t inv_pt_table[N_table_inv_pt];

  static constexpr int ceillog2(int x) { return (x <= 2) ? 1 : 1 + ceillog2((x + 1) / 2); }

  static constexpr int floorlog2(int x) { return (x < 2) ? 0 : 1 + floorlog2(x / 2); }

  static constexpr int pow2(int x) { return x == 0 ? 1 : 2 * pow2(x - 1); }

  template <class data_T, int N>
  static inline float real_val_from_idx(unsigned i) {
    // Treat the index as the top N bits
    static constexpr int NB = ceillog2(N);  // number of address bits for table
    data_T x(0);
    // The MSB of 1 is implicit in the table
    x[x.width - 1] = 1;
    // So we can use the next NB bits for real data
    x(x.width - 2, x.width - NB - 1) = i;
    return (float)x;
  }

  template <class data_T, int N>
  static inline unsigned idx_from_real_val(data_T x) {
    // Slice the top N bits to get an index into the table
    static constexpr int NB = ceillog2(N);  // number of address bits for table
    // Slice the top-1 NB bits of the value
    // the MSB of '1' is implicit, so only slice below that
    ap_uint<NB> y = x(x.width - 2, x.width - NB - 1);
    return (unsigned)y(NB - 1, 0);
  }

  template <class data_T, class table_T, int N>
  static void init_invert_table(table_T table_out[N]) {
    // The template data_T is the data type used to address the table
    for (unsigned i = 0; i < N; i++) {
      float x = real_val_from_idx<data_T, N>(i);
      table_T inv_x = 1 / x;
      table_out[i] = inv_x;
    }
  }

  template <class in_t, class table_t, int N>
  static table_t invert_with_shift(const in_t in, const table_t inv_table[N], bool debug = false) {
    // find the first '1' in the denominator
    int msb = 0;
    for (int b = 0; b < in.width; b++) {
      if (in[b])
        msb = b;
    }
    // shift up the denominator such that the left-most bit (msb) is '1'
    in_t in_shifted = in << (in.width - msb - 1);
    // lookup the inverse of the shifted input
    int idx = idx_from_real_val<in_t, N>(in_shifted);
    table_t inv_in = inv_table[idx];
    // shift the output back
    table_t out = inv_in << (in.width - msb - 1);
    if (debug) {
      std::cout << "           x " << in << ", msb = " << msb << ", shift = " << (in.width - msb) << ", idx = " << idx
                << std::endl;
      std::cout << "     pre 1 / " << in_shifted << " = " << inv_in << "(" << 1 / (float)in_shifted << ")" << std::endl;
      std::cout << "    post 1 / " << in << " = " << out << "(" << 1 / (float)in << ")" << std::endl;
    }
    return out;
  }

  static detaphi_t deltaPhi(Particle a, Particle b) {
    detaphi_t dphi = detaphi_t(a.hwPhi) - detaphi_t(b.hwPhi);
    // phi wrap
    detaphi_t dphi0 = dphi > detaphi_t(l1ct::Scales::INTPHI_PI) ? detaphi_t(l1ct::Scales::INTPHI_TWOPI - dphi) : detaphi_t(dphi);
    detaphi_t dphi1 = dphi < detaphi_t(-l1ct::Scales::INTPHI_PI) ? detaphi_t(l1ct::Scales::INTPHI_TWOPI + dphi) : detaphi_t(dphi);
    detaphi_t dphiw = dphi > detaphi_t(0) ? dphi0 : dphi1;
    return dphiw;
  }

  bool inCone(Particle seed, Particle part) const {
    // scale the particle eta, phi to hardware units
    detaphi_t deta = detaphi_t(seed.hwEta) - detaphi_t(part.hwEta);
    detaphi_t dphi = deltaPhi(seed, part);
    bool ret = deta * deta + dphi * dphi < rCone2;
    //bool ret = r2 < cone2;
    if (_debug) {
      detaphi2_t r2 = detaphi2_t(deta) * detaphi2_t(deta) + detaphi2_t(dphi) * detaphi2_t(dphi);
      std::cout << "  part eta, seed eta: " << part.hwEta << ", " << seed.hwEta << std::endl;
      std::cout << "  part phi, seed phi: " << part.hwPhi << ", " << seed.hwPhi << std::endl;
      std::cout << "  pt, deta, dphi, r2, cone2, lt: " << part.hwPt << ", " << deta << ", " << dphi << ", " << deta * deta + dphi * dphi << ", "
                << rCone2 << ", " << ret << std::endl;
    }
    return ret;
  }

Jet makeJet_HW(const std::vector<Particle>& parts) const {
  // Seed Cone Jet algorithm with ap_fixed types and hardware emulation
  Particle seed = parts.at(0);

  // Fine unless we start using saturation, in which case order matters
  auto sumpt = [](pt_t(a), const Particle& b) { return a + b.hwPt; };

  // Sum the pt
  pt_t pt = std::accumulate(parts.begin(), parts.end(), pt_t(0), sumpt);
  inv_pt_t inv_pt = invert_with_shift<pt_t, inv_pt_t, N_table_inv_pt>(pt, inv_pt_table, false);

  // pt weighted d eta
  std::vector<pt_etaphi_t> pt_deta;
  pt_deta.resize(parts.size());
  std::transform(parts.begin(), parts.end(), pt_deta.begin(), [&seed, &pt](const Particle& part) {
    // In the firmware we calculate the per-particle pt-weighted deta
    return pt_etaphi_t(part.hwPt * detaphi_t(part.hwEta - seed.hwEta));
  });
  // Accumulate the pt-weighted etas. Init to 0, start accumulating at begin()+1 to skip seed
  pt_etaphi_t sum_pt_eta = std::accumulate(pt_deta.begin() + 1, pt_deta.end(), pt_etaphi_t(0));
  etaphi_t eta = seed.hwEta + etaphi_t(sum_pt_eta * inv_pt);

  // pt weighted d phi
  std::vector<pt_etaphi_t> pt_dphi;
  pt_dphi.resize(parts.size());
  std::transform(parts.begin(), parts.end(), pt_dphi.begin(), [&seed, &pt](const Particle& part) {
    // In the firmware we calculate the per-particle pt-weighted dphi
    return pt_etaphi_t(part.hwPt * deltaPhi(part, seed));
  });
  // Accumulate the pt-weighted etas. Init to 0, start accumulating at begin()+1 to skip seed
  pt_etaphi_t sum_pt_phi = std::accumulate(pt_dphi.begin() + 1, pt_dphi.end(), pt_etaphi_t(0));
  etaphi_t phi = seed.hwPhi + etaphi_t(sum_pt_phi * inv_pt);

  Jet jet;
  jet.hwPt = pt;
  jet.hwEta = eta;
  jet.hwPhi = phi;
  jet.constituents = parts;
  /*Jet jet(pt, eta, phi);
  for (auto it = parts.begin(); it != parts.end(); it++) {
    jet.addConstituent(*it);
  }*/

  if (_debug) {
    std::for_each(pt_dphi.begin(), pt_dphi.end(), [](pt_etaphi_t& x) { std::cout << "pt_dphi: " << x << std::endl; });
    std::for_each(pt_deta.begin(), pt_deta.end(), [](pt_etaphi_t& x) { std::cout << "pt_deta: " << x << std::endl; });
    std::cout << " sum_pt_eta: " << sum_pt_eta << ", 1/pt: " << inv_pt << ", sum_pt_eta * 1/pt: " << etaphi_t(sum_pt_eta * inv_pt) << std::endl;
    std::cout << " sum_pt_phi: " << sum_pt_phi << ", 1/pt: " << inv_pt << ", sum_pt_phi * 1/pt: " << etaphi_t(sum_pt_phi * inv_pt) << std::endl;
    std::cout << " uncorr eta: " << seed.hwEta << ", phi: " << seed.hwPhi << std::endl;
    std::cout << "   corr eta: " << eta << ", phi: " << phi << std::endl;
    std::cout << "         pt: " << pt << std::endl;
  }

  return jet;
}

public:
  L1SCJetEmu(bool debug, float coneSize, unsigned nJets) : _debug(debug), _coneSize(coneSize), _nJets(nJets){
    rCone2 = detaphi2_t(coneSize * coneSize / l1ct::Scales::ETAPHI_LSB / l1ct::Scales::ETAPHI_LSB);
    init_invert_table<pt_t, inv_pt_t, N_table_inv_pt>(inv_pt_table);
  }

std::vector<Jet> emulateEvent(std::vector<Particle>& parts) const {
  // The fixed point algorithm emulation
  std::vector<Particle> work;
  work.resize(parts.size());
  std::transform(parts.begin(), parts.end(), work.begin(), [](const Particle& part) { return part; });
  std::sort(work.begin(), work.end(), [](Particle i, Particle j) {
    return (i.hwPt > j.hwPt);
  });

  std::vector<Jet> jets;
  jets.reserve(_nJets);
  while (!work.empty() && jets.size() < _nJets) {
    // Take the first (highest pt) candidate as a seed
    Particle seed = work.at(0);
    // Get the particles within a _coneSize of the seed
    std::vector<Particle> particlesInCone;
    std::copy_if(
        work.begin(), work.end(), std::back_inserter(particlesInCone), [&](const Particle& part) {
          return inCone(seed, part);
        });
    if (_debug) {
      std::cout << "Seed: " << seed.hwPt << ", " << seed.hwEta << ", " << seed.hwPhi << std::endl;
      std::for_each(particlesInCone.begin(), particlesInCone.end(), [&](Particle& part) {
        std::cout << "  Part: " << part.hwPt << ", " << part.hwEta << ", " << part.hwPhi << std::endl;
        inCone(seed, part);
      });
    }
    jets.push_back(makeJet_HW(particlesInCone));
    // remove the clustered particles
    work.erase(
        std::remove_if(work.begin(),
                       work.end(),
                       [&](const Particle& part) { return inCone(seed, part); }),
        work.end());
  }
  return jets;
}

};  // class L1SCJetEmu

#endif
