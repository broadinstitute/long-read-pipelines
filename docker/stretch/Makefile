TAG1 = weisburd/stretch:v0.3.4
TAG2 = weisburd/stretch:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
