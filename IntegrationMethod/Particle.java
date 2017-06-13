package IntegrationMethod;

public class Particle {

    int id;
    double rc;
    double x;
    double y;
    double radius;

    double mass;
    double velX = 0;
    double velY = 0;

    double desiredSpeed;

    boolean isWall = false;
    boolean isOut = false;

    double force;

    public Particle(int id, double x, double y, double velX, double velY, double radius, double mass, double desiredSpeed) {
        this.id = id;
        this.x = x;
        this.y = y;
        this.velX = velX;
        this.velY = velY;
        this.radius = radius;
        this.mass = mass;
        this.desiredSpeed = desiredSpeed;
    }

    public int getId() {
        return id;
    }

    public void setId(int id) {
        this.id = id;
    }

    public double getRc() {
        return rc;
    }

    public void setRc(double rc) {
        this.rc = rc;
    }

    public double getY() {
        return y;
    }

    public void setY(double y) {
        this.y = y;
    }

    public double getX() {
        return x;
    }

    public void setX(double x) {
        this.x = x;
    }

    public double getRadius() {
        return radius;
    }

    public void setRadius(double radius) {
        this.radius = radius;
    }

    public double getMass() {
        return mass;
    }

    public void setMass(double mass) {
        this.mass = mass;
    }

    public double getVelX() {
        return velX;
    }

    public void setVelX(double velX) {
        this.velX = velX;
    }

    public double getVelY() {
        return velY;
    }

    public void setVelY(double velY) {
        this.velY = velY;
    }

    public boolean isWall() {
        return isWall;
    }

    public void setWall(boolean wall) {
        isWall = wall;
    }

    public boolean isOut() {
        return isOut;
    }

    public void setOut(boolean out) {
        isOut = out;
    }

    public double getForce() {
        return force;
    }

    public void setForce(double force) {
        this.force = force;
    }

    public double getDesiredSpeed() {
        return desiredSpeed;
    }

    public void setDesiredSpeed(double desiredSpeed) {
        this.desiredSpeed = desiredSpeed;
    }

    @Override
    public String toString() {
        int vel =  (int)(Math.sqrt(velY*velY)*3);
        return id + "\t" + x + "\t" + y + "\t" + radius + "\t" + mass + "\t" + vel + "\t" +0  + "\t" + 0;
    }
}
