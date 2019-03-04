TAG1 = weisburd/expansionhunter:v2.5.5
TAG2 = weisburd/expansionhunter:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
