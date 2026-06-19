// gpu_mem_probe.h -- source-level GPU peak-memory probe.
//
// Bridge++ alt_Accel uses raw cudaMalloc/cudaFree (no caching allocator /
// memory pool), so the device-wide used = total - free reported by
// cudaMemGetInfo at call time is a faithful high-water mark of real
// allocation (no pool hoarding to inflate it). Sampling in-process at every
// V-cycle catches sub-sampling-period spikes that an external nvidia-smi poll
// (~200 ms cadence) would miss.
//
// Caveats: the figure is device-wide (assumes the GPU is used exclusively by
// this process) and includes the fixed CUDA-context baseline (~300-600 MiB),
// which cancels in a same-device A/B difference. It does NOT instrument the
// allocator: it is a sampler, so it sees the peak only as of the calls made.
#ifndef GPU_MEM_PROBE_H
#define GPU_MEM_PROBE_H

#include <cstddef>
#include <cstdio>
#include <cuda_runtime.h>

namespace gpu_mem_probe {

  inline size_t& peak_used_bytes() { static size_t p = 0; return p; }
  inline size_t& device_total()    { static size_t t = 0; return t; }
  inline long&   n_samples()       { static long   n = 0; return n; }

  // Sample now: update the running max of (total - free).
  inline void sample()
  {
    size_t freeB = 0, totalB = 0;
    if (cudaMemGetInfo(&freeB, &totalB) == cudaSuccess) {
      const size_t used = totalB - freeB;
      if (used > peak_used_bytes()) peak_used_bytes() = used;
      device_total() = totalB;
      ++n_samples();
    }
  }

  // Print the peak (and the device total) so far. tag identifies the call site.
  inline void report(const char* tag)
  {
    const double MiB = 1024.0 * 1024.0;
    std::fprintf(stderr,
                 "[gpu_mem_probe] %s : peak used = %.1f MiB / total %.1f MiB"
                 " (%ld samples)\n",
                 tag,
                 peak_used_bytes() / MiB,
                 device_total()    / MiB,
                 n_samples());
  }

} // namespace gpu_mem_probe

#endif // GPU_MEM_PROBE_H
