using UnityEngine;
using System.Collections;

public class Rigid_Bunny : MonoBehaviour {
    bool Launched = false;
    float Dt = 0.005f;  // 原来是 0.015f
    float Dt_2;
    Vector3 V = new Vector3(0, 0, 0);   // velocity
    Vector3 W = new Vector3(0, 0, 0);   // angular velocity

    const float MassPerVertice = 1;
    float Mass;                                 // mass
    float Mass_INV;
    Vector3 Gravity;
    Matrix4x4 I_ref;                            // reference inertia

    float Linear_decay = 0.999f;                // for velocity decay
    float Angular_decay = 0.98f;
    float Restitution = 0.5f;                   // for collision(mu_n)
    float mu_t = 0.9f;


    Vector3[] Radius;                             // ri
    Vector3 MCenter;                              // 质心
    float EPS = 0.05f;                            // buffer


    // Use this for initialization
    void Start() {
        Dt_2 = Dt / 2;

        Mesh mesh = GetComponent<MeshFilter>().mesh;
        Vector3[] vertices = mesh.vertices;

        // 认为每一个顶点的质量为 1
        float m = MassPerVertice;
        Mass = 0;

        // 因为每一结点的质量相同, 质心就是节点位置的中心
        MCenter = new Vector3(0, 0, 0);
        for (int i = 0; i < vertices.Length; i++) {
            MCenter += vertices[i];
            Mass += m;
            float diag = m * vertices[i].sqrMagnitude;
            I_ref[0, 0] += diag;
            I_ref[1, 1] += diag;
            I_ref[2, 2] += diag;
            I_ref[0, 0] -= m * vertices[i][0] * vertices[i][0];
            I_ref[0, 1] -= m * vertices[i][0] * vertices[i][1];
            I_ref[0, 2] -= m * vertices[i][0] * vertices[i][2];
            I_ref[1, 0] -= m * vertices[i][1] * vertices[i][0];
            I_ref[1, 1] -= m * vertices[i][1] * vertices[i][1];
            I_ref[1, 2] -= m * vertices[i][1] * vertices[i][2];
            I_ref[2, 0] -= m * vertices[i][2] * vertices[i][0];
            I_ref[2, 1] -= m * vertices[i][2] * vertices[i][1];
            I_ref[2, 2] -= m * vertices[i][2] * vertices[i][2];
        }
        I_ref[3, 3] = 1;
        MCenter /= vertices.Length;
        Gravity = new Vector3(0, -Mass * 9.8f, 0);
        Mass_INV = 1.0f / Mass;

        // 计算 ri
        Radius = new Vector3[vertices.Length];
        for (int i = 0; i < vertices.Length; ++i) {
            Radius[i] = vertices[i] - MCenter;
        }
    }

    Matrix4x4 Get_Cross_Matrix(Vector3 a) {
        //Get the cross product matrix of vector a
        Matrix4x4 A = Matrix4x4.zero;
        A[0, 0] = 0;
        A[0, 1] = -a[2];
        A[0, 2] = a[1];
        A[1, 0] = a[2];
        A[1, 1] = 0;
        A[1, 2] = -a[0];
        A[2, 0] = -a[1];
        A[2, 1] = a[0];
        A[2, 2] = 0;
        A[3, 3] = 1;
        return A;
    }

    // In this function, update v and w by the impulse due to the collision with
    // a plane <P, N>
    void Collision_Impulse(Vector3 P, Vector3 N) {
        N = N.normalized;

        // 位置
        Vector3 xHole = transform.position;
        // 旋转矩阵
        Matrix4x4 rHole = Matrix4x4.Rotate(transform.rotation);

        // 获取每一个顶点
        Mesh mesh = GetComponent<MeshFilter>().mesh;
        Vector3[] vertices = mesh.vertices;

        var omegaMatrix = Get_Cross_Matrix(W);
        int hitCount = 0;
        Vector3 hitPos = new Vector3(0, 0, 0);

        // 碰撞检测, 如果有多个点发生碰撞, 则取他们的位置平均
        for (int i = 0; i < vertices.Length; ++i) {
            // Vector4 可以隐式转换为 Vector3（w 被丢弃）
            Vector4 dxi = rHole * new Vector4(Radius[i].x, Radius[i].y, Radius[i].z, 1);
            Vector3 xi = xHole + (Vector3)dxi;

            // 未发生碰撞
            if (Vector3.Dot(xi - P, N) >= EPS) {
                continue;
            }

            // 速度仍然还是向内, 这个时候才需要进行碰撞响应
            Vector3 vi = V + (Vector3)(omegaMatrix * dxi);
            if (Vector3.Dot(vi, N) < 0) {
                ++hitCount;
                hitPos += xi;
            }
        }

        // 没有碰撞发生
        if (hitCount == 0) {
            return;
        }

        // 发生碰撞, 利用平均值求响应
        hitPos /= hitCount;

        // 和上面一样求得 x,v
        Vector3 r = hitPos - xHole;
        Vector4 Rr = new Vector4(r.x, r.y, r.z, 1); // 不需要左乘 rHole 了, 上面乘过了
        Vector3 x = xHole + (Vector3)Rr;
        //float phi = Vector3.Dot(hitPos - P, N);
        Vector3 v = V + (Vector3)(omegaMatrix * Rr);

        // 碰撞响应
        Vector3 vn = Vector3.Dot(v, N) * N;
        Vector3 vt = v - vn;
        float a = Mathf.Max(0, 1 - mu_t * (1 + Restitution) * vn.magnitude / vt.magnitude);
        vn = -Restitution * vn;
        vt = a * vt;
        Vector3 vNew = vt + vn;

        // 计算冲量 j
        Matrix4x4 I_inv = (rHole * I_ref * rHole.transpose).inverse;
        // Matrix4x4 I_inv = I_ref.inverse;
        Matrix4x4 Rrm = Get_Cross_Matrix((Vector3)(Rr));
        Matrix4x4 K = Rrm * I_inv * Rrm;
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j < 4; ++j) {
                K[i, j] = -K[i, j];
            }
            K[i, i] += Mass_INV;
            //K[i, 3] = K[3, i] = 0;
        }
        K[3, 3] = 1;
        Vector4 ft = K.inverse * (vNew - v);

        // 更新 V,W
        V = V + Mass_INV * (Vector3)ft;
        W = W + (Vector3)(I_inv * Rrm * ft);
    }

    // Update is called once per frame
    void Update() {
        // Game Control
        if (Input.GetKey("r")) {
            transform.position = new Vector3(0, 0.6f, 0);
            transform.rotation = new Quaternion();
            Restitution = 0.5f;
            Launched = false;
        }
        if (Input.GetKey("l")) {
            V = new Vector3(5, 2, 0);
            //V = Vector3.zero;
            W = new Vector3(5, 2, 0);
            Launched = true;
        }
        if (!Launched) { return; }

        //////////////////////////////////////////////
        //  Part 0: preparation
        //////////////////////////////////////////////

        // 整体的初始状态
        Quaternion qHole = transform.rotation;
        Vector3 xHole = transform.position;

        // 通用量
        // 旋转矩阵
        Matrix4x4 rHole = Matrix4x4.Rotate(transform.rotation);

        //////////////////////////////////////////////
        // Part I: Update velocities
        //////////////////////////////////////////////

        // 重力影响速度
        V = V + Dt * Mass_INV * Gravity;
        V = Linear_decay * V;
        // 重力不影响角速度
        W = Angular_decay * W;

        //////////////////////////////////////////////
        // Part II: Collision Impulse
        //////////////////////////////////////////////

        Collision_Impulse(new Vector3(0, 0.01f, 0), new Vector3(0, 1, 0));
        Collision_Impulse(new Vector3(2, 0, 0), new Vector3(-1, 0, 0));

        //////////////////////////////////////////////
        // Part III: Update position & orientation
        //////////////////////////////////////////////

        // Update linear status
        xHole = xHole + (Dt * V);
        // Update angular status
        Quaternion q = new Quaternion(Dt_2 * W.z, Dt_2 * W.y, Dt_2 * W.z, 0) * qHole;
        qHole = new Quaternion(q.x + qHole.x, q.y + qHole.y, q.z + qHole.z, q.w + qHole.w);
        // Part IV: Assign to the object
        transform.position = xHole;
        transform.rotation = qHole.normalized;
    }
}
