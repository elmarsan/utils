#ifndef H_TYPES
#define H_TYPES

typedef signed char s8;
typedef unsigned char u8;
typedef short int s16;
typedef unsigned short int u16;
typedef int s32;
typedef unsigned int u32;
typedef signed long long s64;
typedef unsigned long long u64;

#define BIT(x) 1 << (x)
#define KB(x) ((u64)1024 * x)
#define MB(x) ((u64)1024 * KB(x))
#define GB(x) ((u64)1024 * MB(x))

#endif
