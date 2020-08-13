#include <unistd.h> 
#include <stdio.h> 
#include <sys/socket.h> 
#include <stdlib.h> 
#include <netinet/in.h> 
#include <string.h> 
#include <iostream>
#include <array>
#include <netdb.h>
#include <thread>
#include <sys/ioctl.h>
#include <signal.h>
#include <arpa/inet.h> 
#include <string.h> 

class client{
    public:
        int port = 8080;
        int sock = 0;
    int setup_client(){
        int valread; 
        struct sockaddr_in serv_addr;  
        if ((sock = socket(AF_INET, SOCK_STREAM, 0)) < 0) 
            { 
                printf("\n Socket creation error \n"); 
                return 0; 
            }          
        serv_addr.sin_family = AF_INET; 
        serv_addr.sin_port = htons(port);

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
    }

    void run_client(){
        char array[2]; 
        int r = read(sock, &array, 2);
        std::cout << array; 
        close(sock);
    }
    
};

int main(int argc, char *argv[]){
    client c;
    c.setup_client();
    c.run_client();
    return 0;
}