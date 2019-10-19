TAG1 = kgarimella/lr-svim:0.1.0
TAG2 = kgarimella/lr-svim:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
