TAG1 = kgarimella/svim:1.2.0
TAG2 = kgarimella/svim:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
