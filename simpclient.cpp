#include <stdio.h> 
#include <sys/socket.h> 
#include <arpa/inet.h> 
#include <unistd.h> 
#include <string.h> 
#include <stdlib.h>
#define PORT 8080 
   
int main(int argc, char *argv[]) 
{ 
    if (argc != 2) {
        printf("Usage: ./client <num workers>\n");
        exit(2);
    }

    int sock = 0, valread; 
    struct sockaddr_in serv_addr; 
    char *num = argv[1];
    int numprocess = atoi(num);
    int numother = numprocess % 4;
    
    int buffer[1024] = {0}; 
    if ((sock = socket(AF_INET, SOCK_STREAM, 0)) < 0) 
    { 
        printf("\n Socket creation error \n"); 
        return -1; 
    } 
    
    serv_addr.sin_family = AF_INET; 
    serv_addr.sin_port = htons(PORT); 
       
    // Convert IPv4 and IPv6 addresses from text to binary form 
    if(inet_pton(AF_INET, "127.0.0.1", &serv_addr.sin_addr)<=0)  
    { 
        printf("\nInvalid address/ Address not supported \n"); 
        return -1; 
    } 
    
    if (connect(sock, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0) 
    { 
        printf("\nConnection Failed \n"); 
        return -1; 
    } 
    send(sock , &numother , sizeof(numother), 0 ); 
    printf("Nums sent\n"); 
    
    return 0;
} 