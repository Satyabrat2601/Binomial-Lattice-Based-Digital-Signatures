
#ifndef DUET_COMMON_H
#define DUET_COMMON_H

#define DUET_CONCAT_(a, b)  a##b
#define DUET_CONCAT(a, b)   DUET_CONCAT_(a, b)
#define DUET_STR_(x)        #x
#define DUET_STR(x)         DUET_STR_(x)

#endif /* !DUET_COMMON_H */
