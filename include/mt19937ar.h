#ifndef MT19937AR_H
#define MT19937AR_H

/*******************************************************************************************************
 * 🎲 mt19937ar.h
 * -------------------------------------------------------------------------------------------------
 * Mersenne Twister pseudo-random number generator (Matsumoto & Nishimura).
 *******************************************************************************************************/

void init_genrand(unsigned long s);
void init_by_array(unsigned long init_key[], int key_length);

unsigned long genrand_int32(void);
long genrand_int31(void);
double genrand_real1(void);
double genrand_real2(void);
double genrand_real3(void);
double genrand_res53(void);

// Save / restore internal state (for checkpoints)
void get_mt_state(unsigned long *state, int *index);
void set_mt_state(const unsigned long *state, int index);

#endif // MT19937AR_H
