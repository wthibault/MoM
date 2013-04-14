/* 
   osc.cpp
   based on oscpack's SimpleReceive.pp

*/

#include "osc.h"
#include <stdlib.h>
#include <string.h>
#include <stdio.h>

#define PORT 8000 

void
IipPacketListener::handleButton ( const osc::ReceivedMessage& m, const char * name )
{
  float button;
  osc::ReceivedMessage::const_iterator arg = m.ArgumentsBegin();
  if ( arg->IsFloat() )
    button = arg->AsFloat();
  else if ( arg->IsBool() )
    button = arg->AsBool()?1.0:0.0;
  else if ( arg->IsInt32() )
    button = arg->AsInt32();
  pthread_mutex_lock ( &(_params->mutex) );
  if ( _params->value[name] != button ) {
    _params->value[name] = button;
    _params->changed[name] = true;
  }
  pthread_mutex_unlock ( &(_params->mutex) );
}



void
IipPacketListener::handleMidi ( const osc::ReceivedMessage &m )
{
  float  value;
  string addr (m.AddressPattern());
  // ex) /midi/cc40/1
  //     012345678901
  string name = addr.substr(6,4);
  //  cout << "handleMidi " <<  name << endl;
  osc::ReceivedMessage::const_iterator arg = m.ArgumentsBegin();
  if ( arg->IsFloat() )
    value = arg->AsFloat();
  else if ( arg->IsBool() )
    value = arg->AsBool()?1.0:0.0;
  else if ( arg->IsInt32() )
    value = arg->AsInt32();
  pthread_mutex_lock ( &(_params->mutex) );
  if ( _params->value[name] != value ) {
    _params->value[name] = value;
    _params->changed[name] = true;
  }
  pthread_mutex_unlock ( &(_params->mutex) );
  //  cout << "midi at " << name << " with " << value << endl;
}



void
IipPacketListener::handleStringApp ( const osc::ReceivedMessage& m )
{

  string addr (m.AddressPattern());
  osc::ReceivedMessage::const_iterator arg = m.ArgumentsBegin();

  // ex) /StringApp/drawNodes
  //     012345678901234567890
  string name = addr.substr(11,string::npos);

  int ival;
  float fval;
  try {
    ival = arg->AsInt32();
    fval = ival;
  } catch(osc::Exception& e) {
    fval = arg->AsFloat();
  }


  pthread_mutex_lock ( &(_params->mutex) );
  if ( _params->value[name] != fval ) {
    _params->value[name] = fval;
    _params->changed[name] = true;
    std::cout << "handleStringApp " << name << ':' << fval 
	      << "params=" << _params 
	      << std::endl;

  }
  pthread_mutex_unlock ( &(_params->mutex) );

}


void 
IipPacketListener::ProcessMessage( const osc::ReceivedMessage& m, 
                                   const IpEndpointName& remoteEndpoint )
{

  try{
    
    if( strcmp( m.AddressPattern(), "/wii/1/button/A" ) == 0 ){

      handleButton ( m, "A" );
      
    } else if( strcmp( m.AddressPattern(), "/wii/1/button/B" ) == 0 ){

      handleButton ( m, "B" );

    } else if( strcmp( m.AddressPattern(), "/wii/1/button/Home" ) == 0 ){

      handleButton ( m, "Home" );


    } else if( strcmp( m.AddressPattern(), "/wii/1/button/Minus" ) == 0 ){

      handleButton ( m, "Minus" );

    } else if( strcmp( m.AddressPattern(), "/wii/1/button/Plus" ) == 0 ){

      handleButton ( m, "Plus" );

    } else if( strcmp( m.AddressPattern(), "/wii/1/button/Up" ) == 0 ){

      handleButton ( m, "Up" );

    } else if( strcmp( m.AddressPattern(), "/wii/1/button/Down" ) == 0 ){

      handleButton ( m, "Down" );

    } else if( strcmp( m.AddressPattern(), "/wii/1/button/1" ) == 0 ){

      handleButton ( m, "1" );

    } else if( strcmp( m.AddressPattern(), "/wii/1/button/2" ) == 0 ){

      handleButton ( m, "2" );

    } else if( strcmp( m.AddressPattern(), "/wii/1/button/Left" ) == 0 ){

      handleButton ( m, "Left" );

    } else if( strcmp( m.AddressPattern(), "/wii/1/button/Right" ) == 0 ){

      handleButton ( m, "Right" );

    } else if ( string( m.AddressPattern()).find("/midi/cc") != string::npos  ) {

      handleMidi ( m );

    } else if ( string( m.AddressPattern()).find("/StringApp") != string::npos  ) {
      
      handleStringApp ( m );

    } else { }

  } catch ( osc::Exception& e ) {
    std::cout << "error while parsing message: "
              << m.AddressPattern() << ": " << e.what() << "\n";
  }
}

void *
oscThreadFunction ( void * arg )
{
  OscParams *params = static_cast< OscParams* > (arg);
  IipPacketListener listener ( params );
  UdpListeningReceiveSocket s( IpEndpointName( IpEndpointName::ANY_ADDRESS, PORT ),
                               &listener );

    s.RunUntilSigInt();

    return 0;
}

