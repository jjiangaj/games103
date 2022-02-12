using UnityEngine;
using System.Collections;

public class Rigid_Bunny : MonoBehaviour 
{
	bool launched 		= false;
	float dt 			= 0.015f;
	Vector3 v 			= new Vector3(0, 0, 0);	// velocity
	Vector3 w 			= new Vector3(0, 0, 0);	// angular velocity

	float miu_t = 0.5f;

	float mass;									// mass
	Matrix4x4 I_ref;							// reference inertia

	float linear_decay	= 0.999f;				// for velocity decay
	float angular_decay	= 0.98f;				
	float restitution 	= 0.5f;					// for collision
	float g             = 9.8f;


	// Use this for initialization
	void Start () 
	{		
		Mesh mesh = GetComponent<MeshFilter>().mesh;
		Vector3[] vertices = mesh.vertices;

		float m=1;
		mass=0;
		for (int i=0; i<vertices.Length; i++) 
		{
			mass += m;
			float diag=m*vertices[i].sqrMagnitude;
			I_ref[0, 0]+=diag;
			I_ref[1, 1]+=diag;
			I_ref[2, 2]+=diag;
			I_ref[0, 0]-=m*vertices[i][0]*vertices[i][0];
			I_ref[0, 1]-=m*vertices[i][0]*vertices[i][1];
			I_ref[0, 2]-=m*vertices[i][0]*vertices[i][2];
			I_ref[1, 0]-=m*vertices[i][1]*vertices[i][0];
			I_ref[1, 1]-=m*vertices[i][1]*vertices[i][1];
			I_ref[1, 2]-=m*vertices[i][1]*vertices[i][2];
			I_ref[2, 0]-=m*vertices[i][2]*vertices[i][0];
			I_ref[2, 1]-=m*vertices[i][2]*vertices[i][1];
			I_ref[2, 2]-=m*vertices[i][2]*vertices[i][2];
		}
		I_ref [3, 3] = 1;
	}
	
	Matrix4x4 Get_Cross_Matrix(Vector3 a)
	{
		//Get the cross product matrix of vector a
		Matrix4x4 A = Matrix4x4.zero;
		A [0, 0] = 0; 
		A [0, 1] = -a [2]; 
		A [0, 2] = a [1]; 
		A [1, 0] = a [2]; 
		A [1, 1] = 0; 
		A [1, 2] = -a [0]; 
		A [2, 0] = -a [1]; 
		A [2, 1] = a [0]; 
		A [2, 2] = 0; 
		A [3, 3] = 1;
		return A;
	}

	Matrix4x4 Times(Matrix4x4 m, float k)
	{
		Matrix4x4 M = m;
		for(int i = 0; i < 4; i++){
			for(int j = 0; j < 4; j++){
				M [i, j] *= k;
			}
		}
		return M;
	}

	Matrix4x4 MatMinus(Matrix4x4 m, Matrix4x4 n)
	{
		Matrix4x4 A = Matrix4x4.zero;
		for(int i = 0; i < 4; i++){
			for(int j = 0; j < 4; j++){
				A [i, j] =  m [i,j] - n[i, j];
			}
		}
		return A;
	}

	Quaternion QuatAdd(Quaternion s, Quaternion v)
	{
		Quaternion q = new Quaternion(s.x + v.x, s.y + v.y, s.z + v.z, s.w + v.w);
		return q;
	}

	 // In this function, update v and w by the impulse due to the collision with
	 // a plane <P, N>
	 void Collision_Impulse(Vector3 P, Vector3 N)
	 {
	 	Mesh mesh = GetComponent<MeshFilter>().mesh;
	 	Vector3[] vertices = mesh.vertices;
	 	Matrix4x4 R = Matrix4x4.Rotate(transform.rotation);
	 	Vector3 x    = transform.position;
	 	
	 	// Collision detection
	 	Vector3 collidVector3 = new Vector3(0,0,0);
	    Vector3 collidvi            = new Vector3(0, 0, 0);
	 	float count = 0;
	 	
	 	for (int i=0; i<vertices.Length; i++)
	 	{
	 		Vector3 Rri = R.MultiplyPoint(vertices[i]);
	 		if (Vector3.Dot((x + Rri - P), N) < 0)
	        {
		        Vector3 vi = v + Get_Cross_Matrix(w).MultiplyPoint(Rri);
		        if (Vector3.Dot(vi, N) < 0)
		        {
			        collidVector3 += Rri;
			        collidvi  += vi;
			        count++;
		        }
	        }
	 	}
	 	
	 	if (count > 0)
	 	{
	 		// Get average vector
	 		collidVector3 *= (1 / count);
	        collidvi  *= (1 / count);
	 		
	 		Matrix4x4 Rri_cross = Get_Cross_Matrix(collidVector3);

	        Vector3 j = new Vector3(0, 0, 0);
	 
		    Vector3 vni = Vector3.Dot(collidvi , N) * N;
		    Vector3 vti = collidvi  - vni;
		    float a = Mathf.Max(1 - miu_t * (1 + restitution) * vni.magnitude / vti.magnitude, 0); 
		    // 1 - μt(1 + μn)||Vni||/||Vti||
	 			
		    vni *= -restitution;
		    vti *= a;
	 			
		    Matrix4x4 K = MatMinus( Times(Matrix4x4.identity, (1 / mass)) , (Rri_cross * I_ref.inverse * Rri_cross) );
		    j = K.inverse * (vni + vti - collidvi ); // vni+vti = vnew
	            
		    // Update Velocity
	 		v += (1 / mass) * j;
	 		
	 		// Update Angular velocity
	 		w += I_ref.inverse.MultiplyPoint((Rri_cross.MultiplyPoint(j)));
	 
	 		restitution *= 0.5f;
	 	}
	 }


	 // Update is called once per frame
	void Update () 
	{
		//Game Control
		if(Input.GetKey("r"))
		{
			transform.position = new Vector3 (0, 0.6f, 0);
			restitution = 0.5f; //reset restitution because it is decreased in collision calculation
			launched=false;
		}
		if(Input.GetKey("l"))
		{
			v = new Vector3 (5, 2, 0);
			restitution = 0.5f;
			launched=true;
		}

		if (launched != true) return;
		
		// Part I: Update velocities
		Vector3 Fg = new Vector3(0, -mass*g, 0);
		v = v * linear_decay;
		v = v + dt * Fg * (1/mass); //gravity
		w = w * angular_decay; 
		// No torque when there is only gravity

		// Part II: Collision Impulse
		Collision_Impulse(new Vector3(0, 0.01f, 0), new Vector3(0, 1, 0)); 
		Collision_Impulse(new Vector3(2, 0, 0), new Vector3(-1, 0, 0));

		// Part III: Update position & orientation
		//Update linear status
		Vector3 x    = transform.position;
		x = x + dt * v;
		//Update angular status
		Quaternion q = transform.rotation;
		Quaternion dq = new Quaternion(dt*w.x/2, dt*w.y/2, dt*w.z/2, 0.0f);
		q = QuatAdd(q, dq * q);
		q = Quaternion.Normalize(q);
		
        // omg dropbox
		// Part IV: Assign to the object
		transform.position = x;
		transform.rotation = q;
	}
}
