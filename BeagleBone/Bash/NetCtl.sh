#!/bin/bash

if [ "${1,,}" == dhcp ]; then
    cp $SCRIPTSPATH/BeagleBone/Bash/netDHCP /etc/conf.d/net
    /etc/init.d/net.eth0 restart
elif [ "${1,,}" == home ]; then
    cp $SCRIPTSPATH/BeagleBone/Bash/netHome /etc/conf.d/net
    /etc/init.d/net.eth0 restart
elif [ "${1,,}" == ice ]; then
    cp $SCRIPTSPATH/BeagleBone/Bash/netICe /etc/conf.d/net
    /etc/init.d/net.eth0 restart
else
    echo "Usage: "$0" [dhcp | home | ice]"
fi
