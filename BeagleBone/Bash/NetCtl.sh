#!/bin/bash

if [ "${1,,}" == dhcp ]; then
    cp /etc/conf.d/netDHCP /etc/conf.d/net
    /etc/init.d/net.eth0 restart
elif [ "${1,,}" == home ]; then
    cp /etc/conf.d/netHome /etc/conf.d/net
    /etc/init.d/net.eth0 restart
elif [ "${1,,}" == ice ]; then
    cp /etc/conf.d/netICe /etc/conf.d/net
    /etc/init.d/net.eth0 restart
else
    echo "Usage: "$0" [dhcp | home | ice]"
fi
