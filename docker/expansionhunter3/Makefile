TAG1 = weisburd/expansionhunter3:v3.0.0
TAG2 = weisburd/expansionhunter3:latest

all: build push

build:
	docker build -t $(TAG1) -t $(TAG2) .

push:
	docker push $(TAG1)
	docker push $(TAG2)
