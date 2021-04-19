Instrucciones para poder ejecutar los programas:

=========================================================================
=========================================================================

1. Se deben generar dos directorios dentro del directorio
	 donde se contengan los archivos main.f90 y particle_tracking.f90,
	 llamados VecFieldAnim y PTAnim.

2. Se ejecutan los siguientes comandos desde la terminal, en este orden:

:~$ ifort -o VecField.out main.f90 && ./VecField.out
:~$ ifort particle_tracking.f90 -o PTracking.out && ./PTracking.out

3. Para visualizar en gnuplot la animación para el campo de velocidades,
	 se cambia al directorio VecFieldAnim:
	 
	 :~$ cd VecFieldAnim
	 :~$ gnuplot
	 
	 Se declaran las siguientes variables que son para 
	 reducir la cantidad de vectores y su tamaño para una mejor
	 visualización:
	 
	 > set size square; unset key
	 > e = 3; f = 0.25 
	 > l 'anim.gp'

4. Para poder visualizar las partículas, desde el directorio principal:

	:~$ cd PTAnim
	:~$ gnuplot
	
	> set size square
	> set xr [0:10]; set yr [0:10]
	> l 'anim.gp'

=========================================================================
=========================================================================

Para poder visualizar las líneas de corriente:
Una vez que haya terminado de ejecutarse VecField.out, se debe haber 
generado un archivo llamado VecField.dat. 
Se abre entonces gnuplot en el directorio principal,
y se ejecutan los siguientes comandos:

>l 'contour_file.gp'
>p 'contours.dat' w l

=========================================================================
=========================================================================

Para poder visualizar el campo magnético, no es necesario que termine de 
ejecutarse VecField.out, una vez que comienza a ejecutarse, se va a generar
un archivo llamado B0zField.dat.
Se abre entonces gnuplot en el directorio principal, y se ejecuta el 
siguiente comando:

> p 'B0zField.dat' w image

=========================================================================
=========================================================================

Información adicional:
El código calcula diferentes parámetros como son el campo de 
esfuerzos cortantes, el campo de vorticidad, puntos hiperbólicos 
y elípticos definidos através del parámetro Okubo-Weiss (OK), 
la función de corriente (para las líneas de corriente), y la 
magnitud de la velocidad. Exactamente en las líneas 420 y 421 del 
código main.f90 pueden verse las columnas que corresponden a estas
propiedades. Para poder graficarlas, se escribe en gnuplot lo siguiente:

> p 'VecField.dat' u 1:2:# w image

donde # debe reemplazarse por la propiedad que quiera graficarse, 
por ejemplo, si se desea ver el campo de esfuerzos cortantes, 

> p 'VecField.dat' u 1:2:9 w image

=========================================================================
=========================================================================
