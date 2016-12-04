#include <stdio.h>
#include <stdlib.h>
#define SECOND 1
#define THIRD 2
#define W_TYPE SECOND

#if W_TYPE == SECOND
    #define G(x) (if ((x) < 3) return (5) else return(3))
#elif W_TYPE == THIRD
    #define G(x,y) ((x < 3)? 3: 5)
#endif


int main (int argc, char **argv) {
    if(argc != 3){
        puts("Wrong args number");
        exit(0);
    }
    int a;
    int b; 
    b = atoi(argv[1]);
    a = atoi(argv[2]);

    printf("a FUNC b = %d\n", G(a));
    exit(0);
}

