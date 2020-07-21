#include <stdio.h> 
#include <sys/socket.h> 
#include <arpa/inet.h> 
#include <unistd.h> 
#include <string.h> 
#include <stdlib.h>
#include <iostream>
#include <array>
#define PORT 8080 

void run_simpi(int numprocess);  

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
    if(numprocess <= 4){
        run_simpi(numprocess);
        exit(0); 
    }
    int fullWorkstations;
    int remainderCores; 
    if(numprocess % 4 == 0){
        fullWorkstations = numprocess/4;
        remainderCores = 0;
    }
    else{
        fullWorkstations = numprocess/4 + 1;
        remainderCores = numprocess%4;
    }
    char server_ips[2][30] = {"127.0.0.1","10.0.0.1"}; 
    int buffer[1024] = {0}; 
    if ((sock = socket(AF_INET, SOCK_STREAM, 0)) < 0) 
    { 
        printf("\n Socket creation error \n"); 
        return -1; 
    } 
    
    serv_addr.sin_family = AF_INET; 
    serv_addr.sin_port = htons(PORT); 
    /*
    for(int i; i < fullWorkstations; i ++){   
        // Convert IPv4 and IPv6 addresses from text to binary form 

        if(inet_pton(AF_INET, server_ips[i], &serv_addr.sin_addr)<=0)  
        { 
            printf("\nInvalid address/ Address not supported \n"); 
            return -1; 
        } 
        
        if (connect(sock, (struct sockaddr *)&serv_addr, sizeof(serv_addr)) < 0) 
        { 
            printf("\nConnection Failed \n"); 
            return -1; 
        } 
        if(i == fullWorkstations -1){};
        send(sock , &fullWorkstations , sizeof(numother), 0 ); 
        printf("Nums sent\n"); 
    } 
    */  
    return 0;
} 

void run_simpi(int numprocess){
    char progname[100];
    strcpy(progname, "mpi");
    char user[100];
    strcpy(user, "user");
    std::string worker_count_str = std::to_string(numprocess);
    char * args[4] = {progname, user,  const_cast<char*>(worker_count_str.c_str()), NULL};
    if(fork()==0)
        {
        std::cout << progname << "\n";
        std::cout << user << "\n";
        std::cout << numprocess << "\n";
        execv(progname, args);
        }
    return;
}