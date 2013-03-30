/* 
   osc.h
   based on oscpack's SimpleReceive.pp

*/

#ifndef OSC_H
#define OSC_H

#include <iostream>

#include "osc/OscReceivedElements.h"
#include "osc/OscPacketListener.h"
#include "ip/UdpSocket.h"
#include <map>
#include <string>

#include <pthread.h>

using namespace std;

struct OscParams
{
  map<std::string,float> value;
  map<std::string,bool>  changed;
  pthread_mutex_t mutex;
};


class IipPacketListener : public osc::OscPacketListener {
 public:
  IipPacketListener ( OscParams *params )
    { _params = params; };
protected:
  OscParams *_params;

  virtual void handleButton ( const osc::ReceivedMessage& m, const char *name );
  virtual void handleMidi ( const osc::ReceivedMessage &m );
  virtual void handleStringApp ( const osc::ReceivedMessage &m );
  virtual void ProcessMessage( const osc::ReceivedMessage& m, 
                               const IpEndpointName& remoteEndpoint );

};


void *oscThreadFunction ( void * );

#endif
