/*Author: Marco Montecchi
Department of Energy Technologies
ENEA C.R. Casaccia
Roma - Italy

kSEMAW is a workspace for the analysis of
Spectrophotometric (SP), Ellipsometric (ELI) and
Photothermal Deflection Spectroscopy (PDS) measurements


   Copyright (C) 2020  Marco Montecchi

   This program is free software: you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation, either version 3 of the License, or
   (at your option) any later version.

   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program.  If not, see <https://www.gnu.org/licenses/>.
*/
#include <QApplication>
#include "ksemaw.h"

using namespace std;
 
int main(int argc, char *argv[])
{
// avvio k3Dcarspef
    const QByteArray value = qgetenv("USER");
    QString uName=QString::fromLocal8Bit(value);
    QString cmd="/home/"+uName+"/Workspace/qtSource/ksemaw/ksemawf_startup &";
    system(cmd.toStdString().c_str());

 
// avvio GUI ksemaw
    QApplication app(argc, argv);
    ksemaw w;
    w.show();
    return app.exec();
}
