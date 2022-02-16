#!/bin/bash



if [[ $(which docker) && $(docker --version) ]]; then
    echo "Docker installed"
    # command
  else
    echo "Install docker"
    sudo apt-get update
    sudo apt-get install \
    ca-certificates \
    curl \
    gnupg \
    lsb-release
    curl -fsSL https://download.docker.com/linux/ubuntu/gpg | sudo gpg --dearmor -o /usr/share/keyrings/docker-archive-keyring.gpg
    echo \
    "deb [arch=$(dpkg --print-architecture) signed-by=/usr/share/keyrings/docker-archive-keyring.gpg] https://download.docker.com/linux/ubuntu \
    $(lsb_release -cs) stable" | sudo tee /etc/apt/sources.list.d/docker.list > /dev/null
    sudo apt-get update
    sudo apt-get install docker-ce docker-ce-cli containerd.io
    sudo chmod 666 /var/run/docker.sock
    sudo groupadd docker
    sudo usermod -aG docker $USER
fi


# if [[ $(which pgadmin) && $(pgadmin --version) ]]; then
#     echo "pgadmin installed"
#   else
#     sudo curl https://www.pgadmin.org/static/packages_pgadmin_org.pub | sudo apt-key add
#     sudo sh -c 'echo "deb https://ftp.postgresql.org/pub/pgadmin/pgadmin4/apt/$(lsb_release -cs) pgadmin4 main" > /etc/apt/sources.list.d/pgadmin4.list && apt update'
#     sudo apt install pgadmin4
# fi


if [[ $(which docker-compose) && $(docker-compose --version) ]]; then
    echo "pgadmin installed"
  else
    sudo curl -L "https://github.com/docker/compose/releases/download/1.29.2/docker-compose-$(uname -s)-$(uname -m)" -o /usr/local/bin/docker-compose
    sudo chmod +x /usr/local/bin/docker-compose
    docker-compose --version
fi