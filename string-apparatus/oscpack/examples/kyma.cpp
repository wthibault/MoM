/* 
    Simple example of sending an OSC message using oscpack.
*/

#include "osc/OscOutboundPacketStream.h"
#include "ip/UdpSocket.h"


#define ADDRESS "192.168.2.3"
#define PORT 8000

#define OUTPUT_BUFFER_SIZE 1024

int main(int argc, char* argv[])
{
    UdpTransmitSocket transmitSocket( IpEndpointName( ADDRESS, PORT ) );
    
    char buffer[OUTPUT_BUFFER_SIZE];
    osc::OutboundPacketStream p( buffer, OUTPUT_BUFFER_SIZE );
    
    p << 
        << osc::BeginMessage( "/osc/respond_to" ) 
            << 7000 
      << osc::EndMessage
      ;
    
    transmitSocket.Send( p.Data(), p.Size() );
}

