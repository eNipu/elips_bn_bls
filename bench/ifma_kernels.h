/* Kernels measured by bench/ifma_probe.c. See that file for what and why. */
#ifndef ELIPS_IFMA_KERNELS_H
#define ELIPS_IFMA_KERNELS_H
#include <stdint.h>
uint64_t part_convert(const uint64_t in[6]);
uint64_t part_loop(const uint64_t a52[8], const uint64_t b52[8]);
uint64_t part_loop_nodep(const uint64_t a52[8], const uint64_t b52[8]);
void     part_to52(uint64_t o[8], const uint64_t in[6], unsigned s);
uint64_t part_loop_x2(const uint64_t a[8], const uint64_t b[8]);
uint64_t part_loop_x4(const uint64_t a[8], const uint64_t b[8]);
uint64_t part_full_x4(const uint64_t a[8], const uint64_t b[8]);
void part_from52(uint64_t o[6], const uint64_t in[8]);
#endif
