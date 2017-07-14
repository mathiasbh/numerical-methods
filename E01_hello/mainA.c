#include <stdio.h>
#include <stdlib.h>

// the program is run on my own machine where my user is "mathias"
int main() {
  char *user = getenv("USER"); // looks in environment for USER
  printf("hello, %s\n", user); 
  return 0;
}
