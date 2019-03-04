TAG1 = weisburd/repeat-finder:v1
TAG2 = weisburd/repeat-finder:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
