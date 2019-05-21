#include <stdio.h>
#include <unistd.h> // notice this! you need it!

int main(int argc, char* argv[]){
    printf("Hello,");
    int x= argv[1] - '0';
    sleep(x); // format is sleep(x); where x is # of seconds.
    printf("World");
    return 0;
}
