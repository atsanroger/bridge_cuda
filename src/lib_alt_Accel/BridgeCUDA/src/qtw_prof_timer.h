/*!
      @file    qtw_prof_timer.h
      @brief   Lightweight env-gated per-kernel GPU timer for the QTW DWF
               kernels. nsys/ncu can't time kernels under WSL2 (no CUPTI
               kernel activity), so we use cudaEvent timing around each launch.
               Active only when the env var QTW_PROFILE is set; otherwise the
               launch runs untouched (no events, no sync added beyond the
               wrapper's existing cudaDeviceSynchronize()).
      @author  Wei-Lun Chen (wlchen)
*/
#ifndef QTW_PROF_TIMER_INCLUDED
#define QTW_PROF_TIMER_INCLUDED

#include <cuda_runtime.h>
#include <cstdio>
#include <cstdlib>
#include <map>
#include <string>

// One registry per program (inline function-local static => single instance
// across all translation units that include this header). Its destructor dumps
// the accumulated per-kernel totals at program exit when QTW_PROFILE is set.
struct QtwProfReg {
    std::map<std::string, double> ms;
    std::map<std::string, long>   n;
    ~QtwProfReg() {
        if (!getenv("QTW_PROFILE")) return;
        double tot = 0;
        for (auto& kv : ms) tot += kv.second;
        fprintf(stderr, "\n=== QTW per-kernel GPU time (cudaEvent) ===\n");
        fprintf(stderr, "%12s %7s %8s   kernel\n", "ms", "%", "calls");
        for (auto& kv : ms) {
            fprintf(stderr, "%12.3f %6.1f%% %8ld   %s\n",
                    kv.second, tot > 0 ? 100.0 * kv.second / tot : 0.0,
                    n[kv.first], kv.first.c_str());
        }
        fprintf(stderr, "%12.3f %6.1f%%            TOTAL\n", tot, 100.0);
    }
};

inline QtwProfReg& qtw_prof() { static QtwProfReg r; return r; }
inline bool qtw_prof_on() { static bool on = (getenv("QTW_PROFILE") != nullptr); return on; }
inline void qtw_prof_add(const char* name, float ms) {
    qtw_prof().ms[name] += ms;
    qtw_prof().n[name]  += 1;
}

// Bracket a kernel launch with cudaEvent timing (the launch itself has too
// many top-level commas to pass as a macro argument). Use as:
//     QTW_PROF_BEGIN();
//     my_kernel<<<g,b>>>(...);
//     QTW_PROF_END("my_kernel");
//     CHECK(cudaDeviceSynchronize());
// One pair per function scope (declares _qa/_qb/_qp).
#define QTW_PROF_BEGIN()                                              \
    cudaEvent_t _qa, _qb; bool _qp = qtw_prof_on();                  \
    if (_qp) { cudaEventCreate(&_qa); cudaEventCreate(&_qb);          \
               cudaEventRecord(_qa); }

#define QTW_PROF_END(name)                                           \
    if (_qp) { cudaEventRecord(_qb); cudaEventSynchronize(_qb);       \
        float _qms = 0.f; cudaEventElapsedTime(&_qms, _qa, _qb);      \
        qtw_prof_add(name, _qms);                                     \
        cudaEventDestroy(_qa); cudaEventDestroy(_qb); }

#endif
