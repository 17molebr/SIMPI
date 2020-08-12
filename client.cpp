#include <stdio.h> 
#include <sys/socket.h> 
#include <arpa/inet.h> 
#include <unistd.h> 
#include <string.h> 
#include <stdlib.h>
#include <iostream>
#include <array>
#include <signal.h>
#define PORT 8080 

int main(){
   int sock = 0, valread; 
   struct sockaddr_in serv_addr;  
   if ((sock = socket(AF_INET, SOCK_STREAM, 0)) < 0) 
        { 
            printf("\n Socket creation error \n"); 
            return 0; 
        } 
    
        
    serv_addr.sin_family = AF_INET; 
    serv_addr.sin_port = htons(PORT);

    if(inet_pton(AF_INET, "192.168.168.13", &serv_addr.sin_addr)<=0)  
    { 
        printf("\nInvalid address/ Address not supported \n"); 
        return 0; 
    }
    if (connect(sock, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0) 
    { 
        printf("\nConnection Failed \n"); 
        return 0; 
    } 
    char array[2]; 
    int r = read(sock, &array, 2);
    std::cout << array; 
    close(sock);
    return 0;
}
