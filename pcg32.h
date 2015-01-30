/*
 * Tiny self-contained version of the PCG Random Number Generation for C++
 * put together from pieces of the much larger C codebase by Wenzel Jakob.
 *
 * The PCG random number generator was developed by Melissa O'Neill <oneill@pcg-random.org>
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 *
 * For additional information about the PCG random number generation scheme,
 * including its license and other licensing options, visit
 *
 *     http://www.pcg-random.org
 */

#ifndef PCG32_H_INCLUDED
#define PCG32_H_INCLUDED 1

#include <inttypes.h>
#include <cmath>

/// PCG32 Pseudorandom number generator
struct pcg32 {
	/// Initialize the pseudorandom number generator with default seed
	inline pcg32(uint64_t state = 0x853c49e6748fea9bULL,
	             uint64_t inc = 0xda3e39cb94b95bdbULL) : state(state), inc(inc) { }

	/**
	 * \brief Seed the pseudorandom number generator
	 *
	 * Specified in two parts: a state initializer and a sequence selection
	 * constant (a.k.a. stream id)
	 */
	inline void seed(uint64_t initstate, uint64_t initseq = 1) {
		state = 0U;
		inc = (initseq << 1u) | 1u;
		nextUInt();
		state += initstate;
		nextUInt();
	}

	/// Generate a uniformly distributed unsigned 32-bit random number
	inline uint32_t nextUInt() {
		uint64_t oldstate = state;
		state = oldstate * 6364136223846793005ULL + inc;
		uint32_t xorshifted = (uint32_t) (((oldstate >> 18u) ^ oldstate) >> 27u);
		uint32_t rot = (uint32_t) (oldstate >> 59u);
		return (xorshifted >> rot) | (xorshifted << ((-rot) & 31));
	}

	/// Generate a uniformly distributed number, r, where 0 <= r < bound
	inline uint32_t nextUInt(uint32_t bound) {
		// To avoid bias, we need to make the range of the RNG a multiple of
		// bound, which we do by dropping output less than a threshold.
		// A naive scheme to calculate the threshold would be to do
		//
		//     uint32_t threshold = 0x100000000ull % bound;
		//
		// but 64-bit div/mod is slower than 32-bit div/mod (especially on
		// 32-bit platforms).  In essence, we do
		//
		//     uint32_t threshold = (0x100000000ull-bound) % bound;
		//
		// because this version will calculate the same modulus, but the LHS
		// value is less than 2^32.

		uint32_t threshold = -bound % bound;

		// Uniformity guarantees that this loop will terminate.  In practice, it
		// should usually terminate quickly; on average (assuming all bounds are
		// equally likely), 82.25% of the time, we can expect it to require just
		// one iteration.  In the worst case, someone passes a bound of 2^31 + 1
		// (i.e., 2147483649), which invalidates almost 50% of the range.  In 
		// practice, bounds are typically small and only a tiny amount of the range
		// is eliminated.
		for (;;) {
			uint32_t r = nextUInt();
			if (r >= threshold)
				return r % bound;
		}
	}

	/// Generate a single precision floating point value on the interval [0, 1)
	inline float nextFloat() {
		return std::ldexp((float) nextUInt(), -32);
	}

	/// Generate a double precision floating point value on the interval [0, 1)
	inline double nextDouble() {
		return std::ldexp((double) nextUInt(), -32);
	}

	/// Multi-step advance function (jump-ahead, jump-back)
	inline void advance(int64_t delta) {
		state = pcg_advance_lcg_64(state, (uint64_t) delta,
			6364136223846793005ULL, inc);
	}

	/**
	 * Draw uniformly distributed permutation and permute the
	 * given STL container (Knuth, TAoCP Vol. 2 (3rd 3d), Section 3.4.2)
	 */
	template <typename Iterator> inline void shuffle(Iterator begin, Iterator end) {
		for (Iterator it = end - 1; it > begin; --it)
			std::iter_swap(it, begin + nextUInt((uint32_t) (it - begin + 1)));
	}

	uint64_t state;  // RNG state.  All values are possible.
	uint64_t inc;    // Controls which RNG sequence (stream) is selected. Must *always* be odd.
private:

	/**
	 * Multi-step advance functions (jump-ahead, jump-back)
	 *
	 * The method used here is based on Brown, "Random Number Generation
	 * with Arbitrary Stride,", Transactions of the American Nuclear
	 * Society (Nov. 1994).  The algorithm is very similar to fast
	 * exponentiation.
	 *
	 * Even though delta is an unsigned integer, we can pass a
	 * signed integer to go backwards, it just goes "the long way round".
	 */
	static uint64_t pcg_advance_lcg_64(uint64_t state, uint64_t delta, uint64_t cur_mult, uint64_t cur_plus) {
		uint64_t acc_mult = 1u;
		uint64_t acc_plus = 0u;
		while (delta > 0) {
			if (delta & 1) {
				acc_mult *= cur_mult;
				acc_plus = acc_plus * cur_mult + cur_plus;
			}
			cur_plus = (cur_mult + 1) * cur_plus;
			cur_mult *= cur_mult;
			delta /= 2;
		}
		return acc_mult * state + acc_plus;
	}
};

#endif // PCG32_H_INCLUDED
